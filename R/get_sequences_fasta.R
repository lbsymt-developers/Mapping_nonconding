load("data/cred_promoters_GRanges.rda")

library(GenomicRanges)
apoe <- credpromoter[credpromoter$gene_name %in% "APOE"]
apoe <- apoe[!duplicated(apoe$rsid)]
seqlevels(apoe) <- paste0("chr", seqlevels(apoe))
library(BSgenome.Hsapiens.UCSC.hg38)
seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, apoe)
# Aqui ponerle el nombre de los SNPs de donde provienen
names(seq) <- paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, "APOE_promoters.fasta")

