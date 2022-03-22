promoter <- read.table("data/Gencode19_promoter.bed")
promoterranges <- GRanges(promoter[,1], IRanges(promoter[,2], promoter[,3]), gene = promoter[,4])
# Overlap credible SNPs with promoter regions.
olap <- findOverlaps(credranges, promoterranges)
credpromoter <- credranges[queryHits(olap)]
mcols(credpromoter) <- cbind(mcols(credpromoter), mcols(promoterranges[subjectHits(olap)]))

gen_anno <- vroom::vroom("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit")
gen_anno <- gen_anno[!is.na(gen_anno$`Ensembl gene ID`), ]
save(gen_anno, file = "data/Gen_Anno.rda")

