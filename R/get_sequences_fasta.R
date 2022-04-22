load("data/promoter_ranges.rda")


apoe <- promoterranges[promoterranges$gene_name %in% "APOE"]
seqlevels(apoe) <- paste0("chr", seqlevels(apoe))
library(BSgenome.Hsapiens.UCSC.hg38)
seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, apoe)
# Aqui ponerle el nombre de los SNPs de donde provienen
names(seq) <- paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, "my.fasta")

