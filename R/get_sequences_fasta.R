load("data/cred_promoters_GRanges.rda")

library(BSgenome.Hsapiens.UCSC.hg38)
seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, olaps)
# Aqui ponerle el nombre de los SNPs de donde provienen
names(seq) <- paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, "my.fasta")
