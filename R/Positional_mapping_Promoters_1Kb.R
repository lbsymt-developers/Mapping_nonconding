library(GenomicRanges)

promoter <- read.table("data/gencodev39_promoters1000pb.bed.txt")
ens <- stringr::str_split(promoter$V4, "_up")
ens <- sapply(ens, function(x){
  x[[1]]
})

promoter$V4 <- ens

genes <- stringr::str_split(promoter$V1, "chr")
genes <- sapply(genes, function(x){
  x[[2]]
})

genes <- stringr::str_split(genes, "_")
genes <- sapply(genes, function(x){
  x[[1]]
})

promoter$V1 <- genes

# AHORA CONVERTIMOS A ENSG

library(AnnotationHub)
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v104")
edb <- ah[["AH95744"]]
ids <- ens
txs <- transcripts(edb, columns = c("tx_id", "tx_biotype",
                                    "tx_id_version", "gc_content",
                                    "gene_name", "gene_id"))
idx <- match(ids, txs$tx_id_version)
idx <- na.omit(idx)
txs_promoter <- txs[idx, ]
# save(txs_promoter, file = "txs_promoters_annotation.rda")
df_txs <- data.frame(txs_promoter)

colnames(promoter)[4] <- "tx_id_version"
promoters_txs <- dplyr::left_join(promoter, df_txs, "tx_id_version")

promoterranges <- GRanges(promoters_txs[,1],
                          IRanges(promoters_txs[,2], promoters_txs[,3]),
                          Tx_id_version = promoters_txs[,4], tx_biotype = promoters_txs$tx_biotype,
                          gene_name = promoters_txs$gene_name,
                          gene_id = promoters_txs$gene_id,
                          tx_id = promoters_txs$tx_id)

# save(promoterranges, file = "data/promoter_ranges.rda")


load("data/AD_credibleSNP.rda")
# Overlap credible SNPs with promoter regions.
olap <- findOverlaps(credranges, promoterranges)
credpromoter <- credranges[queryHits(olap)]
mcols(credpromoter) <- cbind(mcols(credpromoter), mcols(promoterranges[subjectHits(olap)]))
# save(credpromoter, file = "data/cred_promoters_GRanges.rda")

load("data/cred_promoters_GRanges.rda")
df <- data.frame(credpromoter)
df <- na.omit(df)

table(df$gene_name)
summary(df$gene_name)

# Se tienen 56 genes distintos mapeados

# "https://www.genenames.org/download/custom/"
# gen_anno <- vroom::vroom("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit")
# gen_anno <- gen_anno[!is.na(gen_anno$`Ensembl gene ID`), ]
# save(gen_anno, file = "data/Gen_Anno.rda")

genes <- df$gene_name
gen_sel <- genes[genes != ""]
gen_sel <- gen_sel[!duplicated(gen_sel)]

# dir.create("resultados")
saveRDS(gen_sel, file = "resultados/genes_snps_promotores.rds")
