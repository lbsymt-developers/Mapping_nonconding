load("data/promoter_ranges.rda")
load("resultados/cred_promoters_GRanges.rda")

library(GenomicRanges)
apoe <- credpromoter[credpromoter$gene_name %in% "APOE"]
apoe <- apoe[!duplicated(apoe$rsid)]
seqlevels(apoe) <- paste0("chr", seqlevels(apoe))
library(BSgenome.Hsapiens.UCSC.hg38)
promoters_apoe <- GRanges(seqnames = "chr19", IRanges(start= 44904790,
                                                      end =44905790,
                                                      name = "APOE"))
seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, promoters_apoe)
# Aqui ponerle el nombre de los SNPs de donde provienen
names(seq) <- paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, "APOE_promoters.fasta")

subseq(seq, start = 558, end = 573)

posi <- c(789, 581, 517)
df$pos <- posi
