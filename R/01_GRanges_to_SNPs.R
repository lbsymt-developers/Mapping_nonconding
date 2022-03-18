library(GenomicRanges)
options(stringsAsFactors = F)
credSNP <- vroom::vroom("data/Supplementary_Table_8_Jansen.csv")
credSNP <- credSNP[credSNP$`Credible Causal`=="Yes",]

credranges <- GRanges(credSNP$Chr, IRanges(credSNP$bp, credSNP$bp), rsid = credSNP$SNP, P=credSNP$P)
save(credranges, file ="data/AD_credibleSNP.rda")

# "http://genome.ucsc.edu/cgi-bin/hgTables"
