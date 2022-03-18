library(GenomicRanges)
options(stringsAsFactors = F)
credSNP <- vroom::vroom("datos_newAnalysis/AD_complete.csv")
# credSNP <- credSNP[credSNP$`Credible Causal`=="Yes",]

credranges <- GRanges(credSNP$seq_region_name, IRanges(credSNP$start, credSNP$start),
                      rsid = credSNP$SNP, Context=credSNP$CONTEXT)


# Load promoter and exonic region and generate a GRange object.
exon <- read.table("data/Gencode19_exon.bed")
exonranges <- GRanges(exon[,1],IRanges(exon[,2], exon[,3]), gene = exon[,4])
promoter <- read.table("data/Gencode19_promoter.bed")
promoterranges <- GRanges(promoter[,1], IRanges(promoter[,2], promoter[,3]), gene = promoter[,4])

# Overlap credible SNPs with exonic regions.
olap <- findOverlaps(credranges, exonranges)
credexon <- credranges[queryHits(olap)]
mcols(credexon) <- cbind(mcols(credexon), mcols(exonranges[subjectHits(olap)]))

# Overlap credible SNPs with promoter regions.
olap <- findOverlaps(credranges, promoterranges)
credpromoter <- credranges[queryHits(olap)]
mcols(credpromoter) <- cbind(mcols(credpromoter), mcols(promoterranges[subjectHits(olap)]))

# Link SNPs to their putative target genes using chromatin interactions.
# Load Hi-C dataset and generate a GRange object.
hic <- read.table("data/Promoter-anchored_chromatin_loops.bed", skip=1)
colnames(hic) <- c("chr", "TSS_start", "TSS_end", "Enhancer_start", "Enhancer_end")
hicranges <- GRanges(hic$chr, IRanges(hic$TSS_start, hic$TSS_end), enhancer=hic$Enhancer_start)
olap <- findOverlaps(hicranges, promoterranges)
hicpromoter <- hicranges[queryHits(olap)]
mcols(hicpromoter) <- cbind(mcols(hicpromoter), mcols(promoterranges[subjectHits(olap)]))
hicenhancer <- GRanges(seqnames(hicpromoter), IRanges(hicpromoter$enhancer, hicpromoter$enhancer+10000),
                       gene=hicpromoter$gene)

# Overlap credible SNPs with Hi-C GRange object.
olap <- findOverlaps(credranges, hicenhancer)
credhic <- credranges[queryHits(olap)]
mcols(credhic) <- cbind(mcols(credhic), mcols(hicenhancer[subjectHits(olap)]))

# Compile AD candidate genes defined by positional mapping and chromatin interaction profiles.
### The resulting candidate genes for AD:
ADgenes <- Reduce(union, list(credhic$gene, credexon$gene, credpromoter$gene))

SNPs <- c(credhic$rsid, credexon$rsid, credpromoter$rsid)
SNPs <- unique(SNPs)
saveRDS(SNPs, "datos_newAnalysis/SNPs.rds")

### to convert Ensembl Gene ID to HGNC symbol
load("data/geneAnno2.rda")
ADhgnc <- geneAnno1[match(ADgenes, geneAnno1$ensembl_gene_id), "hgnc_symbol"]
ADhgnc <- ADhgnc[ADhgnc!=""]
# save(ADgenes, ADhgnc, file="datos_newAnalysis/ADgenes.rda")
# write.table(ADhgnc, file="datos_newAnalysis/ADgenes.txt", row.names=F, col.names=F, quote=F, sep="\t")

library(reshape)
library(ggplot2)
library(GenomicRanges)
library(biomaRt)
library("WGCNA")

# Process expression and meta data
datExpr <- read.csv("data/gene_array_matrix_csv/expression_matrix.csv", header = FALSE)
datExpr <- datExpr[,-1]
datMeta <- read.csv("data/gene_array_matrix_csv/columns_metadata.csv")
datProbes <- read.csv("data/gene_array_matrix_csv/rows_metadata.csv")
datExpr <- datExpr[datProbes$ensembl_gene_id!="",]
datProbes <- datProbes[datProbes$ensembl_gene_id!="",]
datExpr.cr <- collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id,
                           rowID= rownames(datExpr))
datExpr <- datExpr.cr$datETcollapsed
gename <- data.frame(datExpr.cr$group2row)
rownames(datExpr) <- gename$group

# Specify developmental stages.
datMeta$Unit <- "Postnatal"
idx <- grep("pcw", datMeta$age)
datMeta$Unit[idx] <- "Prenatal"
idx <- grep("yrs", datMeta$age)
datMeta$Unit[idx] <- "Postnatal"
datMeta$Unit <- factor(datMeta$Unit, levels=c("Prenatal", "Postnatal"))

# Select Cortical regions
datMeta$Region <- "SubCTX"
r <- c("A1C", "STC", "ITC", "TCx", "OFC", "DFC",
       "VFC", "MFC", "M1C", "S1C", "IPC", "M1C-S1C",
       "PCx", "V1C", "Ocx")
datMeta$Region[datMeta$structure_acronym %in% r] <- "CTX"
datExpr <- datExpr[,which(datMeta$Region=="CTX")]
datMeta <- datMeta[which(datMeta$Region=="CTX"),]
save(datExpr, datMeta, file="datos_newAnalysis/devExpr.rda")

# Extract developmental expression profiles of AD risk genes.
load("datos_newAnalysis/ADgenes.rda")
exprdat <- apply(datExpr[match(ADgenes, rownames(datExpr)),], 2, mean, na.rm=T)
dat <- data.frame(Region=datMeta$Region, Unit=datMeta$Unit, Expr=exprdat)

prueba <- datExpr[match(ADgenes, rownames(datExpr)),]
prueba <- na.omit(prueba)
ADhgnc <- geneAnno1[match(rownames(prueba), geneAnno1$ensembl_gene_id), "hgnc_symbol"]
ADhgnc <- ADhgnc[ADhgnc!=""]

library(RColorBrewer)
library(gplots)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("cyan","red")[dat$Unit]

# Plot the heatmap
tiff("High_var_genes.heatmap.tiff", height = 9, width = 12,
     units = 'in', res=300, compression = "none")
heatmap.2(prueba,col=rev(morecols(50)),trace="none",
          main="Genes de riesgo AD entre las muestras",
          ColSideColors=col.cell,scale="row")
dev.off()


# Compare prenatal versus postnatal expression levels of AD risk genes.
pdf(file="developmental_expression.pdf")
ggplot(dat,aes(x=Unit, y=Expr, fill=Unit, alpha=Unit)) +
  ylab("Normalized expression") + geom_boxplot(outlier.size = NA) +
  ggtitle("Expresi??n en cerebro: regi??n cortical") + xlab("") +
  scale_alpha_manual(values=c(0.2, 1)) + theme_classic() +
  theme(legend.position="na", text = element_text(size = 20))
dev.off()



load("datos_newAnalysis/ADgenes.rda")
load("data/geneAnno2.rda")
targetname <- "AD"
targetgene <- ADhgnc
cellexp <- read.table("data/DER-20_Single_cell_expression_processed_TPM_backup.tsv",header=T,fill=T)
cellexp[1121,1] <- cellexp[1120,1]
cellexp <- cellexp[-1120,]
rownames(cellexp) <- cellexp[,1]
cellexp <- cellexp[,-1]
datExpr <- scale(cellexp,center=T, scale=F)
datExpr <- datExpr[,789:ncol(datExpr)]

# Extract cellular expression profiles of AD risk genes
exprdat <- apply(datExpr[match(targetgene, rownames(datExpr)),],2,mean,na.rm=T)
dat <- data.frame(Group=targetname, cell=names(exprdat), Expr=exprdat)
dat$celltype <- unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))
dat <- dat[-grep("Ex|In",dat$celltype),]
dat$celltype <- gsub("Dev","Fetal",dat$celltype)
dat$celltype <- factor(dat$celltype, levels=c("Neurons","Astrocytes","Microglia","Endothelial",
                                              "Oligodendrocytes","OPC","Fetal"))
pdf(file ="singlecell_expression_ADgenes.pdf")
ggplot(dat,aes(x=celltype, y=Expr, fill=celltype)) +
  ylab("Normalized expression") + xlab("") +
  geom_violin() + theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust=1),
        text = element_text(size=15)) +
  theme(legend.position="none") +
  ggtitle(paste0("Cellular expression profiles of AD risk genes"))
dev.off()




credSNP <- vroom::vroom("datos_newAnalysis/AD_complete.csv")
SNPs <- readRDS("datos_newAnalysis/SNPs.rds")

seleccionados <- credSNP[credSNP$SNP %in% SNPs, ]
graf <- table(seleccionados$seq_region_name)
graf <- data.frame(chr = names(graf), value = graf)
graf$chr <- as.numeric(graf$chr)
graf$chr <- as.factor(graf$chr)

ggplot(graf, aes(x = chr, y = value.Freq)) +
  geom_col(fill = "darkred", color = "black") + theme_bw() +
  ylab("Number of SNPs") + xlab("Chr") +
  theme(text = element_text(size = 20))

