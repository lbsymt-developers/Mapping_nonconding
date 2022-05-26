# Link SNPs to their putative target genes using chromatin interactions.
# Load Hi-C dataset and generate a GRange object.
library(GenomicRanges)
hic <- read.table("data/Promoter-anchored_chromatin_loops.bed", skip=1)
load("data/promoter_ranges.rda")
colnames(hic) <- c("chr", "TSS_start", "TSS_end", "Enhancer_start", "Enhancer_end")
hicranges <- GRanges(hic$chr, IRanges(hic$TSS_start, hic$TSS_end), enhancer=hic$Enhancer_start)
olap <- GenomicRanges::findOverlaps(hicranges, promoterranges)
hicpromoter <- hicranges[queryHits(olap)]
mcols(hicpromoter) <- cbind(mcols(hicpromoter), mcols(promoterranges[subjectHits(olap)]))
hicenhancer <- GRanges(seqnames(hicpromoter), IRanges(hicpromoter$enhancer, hicpromoter$enhancer+10000),
                       gene=hicpromoter$gene_name,
                       TSS = ranges(hicpromoter),
                       Enhancer = IRanges(hicpromoter$enhancer, hicpromoter$enhancer+10000))

# saveRDS(hicenhancer, "resultados/hicenhancer_granges.rda")

load("data/AD_credibleSNP.rda")
# Overlap credible SNPs with Hi-C GRange object.
olap <- findOverlaps(credranges, hicenhancer)
credhic <- credranges[queryHits(olap)]
mcols(credhic) <- cbind(mcols(credhic), mcols(hicenhancer[subjectHits(olap)]))
saveRDS(credhic, "resultados/snps_mapping_HiC_granges.rds")


df <- data.frame(credhic)
df <- na.omit(df)
df <- df[!duplicated(df$rsid), ]

table(df$gene)
summary(df$gene_name)

# Se tienen 56 genes distintos mapeados
# Hay 403 SNPs no codificantes en 60 genes y 25 pseudogenes y lncRNA

genes <- df$gene
gen_sel <- genes[genes != ""]
gen_sel <- gen_sel[!duplicated(gen_sel)]

# dir.create("resultados")
saveRDS(gen_sel, file = "resultados/genes_snps_HiC.rds")
