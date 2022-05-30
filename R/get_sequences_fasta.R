load("data/promoter_ranges.rda")

library(GenomicRanges)
apoe <- promoterranges[promoterranges$gene_name %in% "APOE"]
# apoe <- apoe[!duplicated(apoe$rsid)]
seqlevels(apoe) <- paste0("chr", seqlevels(apoe))
library(BSgenome.Hsapiens.UCSC.hg38)
promoters_apoe <- GRanges(seqnames = "chr19",
                          IRanges(start = c(44904790, 44904790, 44904790),
                                  end = c(44905790, 44905790, 44905790),
                                  names = c("APOE", "APOE", "APOE")))
seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, apoe)
# Aqui ponerle el nombre de los SNPs de donde provienen
names(seq) <- paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, "APOE_promoters.fasta")

