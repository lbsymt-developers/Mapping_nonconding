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

promoterranges <- GRanges(promoter[,1], IRanges(promoter[,2], promoter[,3]), Tx = promoter[,4])

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

# Overlap credible SNPs with promoter regions.
olap <- findOverlaps(credranges, promoterranges)
credpromoter <- credranges[queryHits(olap)]
mcols(credpromoter) <- cbind(mcols(credpromoter), mcols(promoterranges[subjectHits(olap)]))

# "https://www.genenames.org/download/custom/"
# gen_anno <- vroom::vroom("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit")
# gen_anno <- gen_anno[!is.na(gen_anno$`Ensembl gene ID`), ]
# save(gen_anno, file = "data/Gen_Anno.rda")

