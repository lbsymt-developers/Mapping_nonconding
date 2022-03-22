library(GenomicRanges)
options(stringsAsFactors = F)

gwascatalog_AD <- seekerBio::seeker_gwas("Alzheimer")
gwascatalog_AD <- gwascatalog_AD[gwascatalog_AD$P.VALUE <= 5E-8, ]
coding <- c("missense_variant", "stop_gained", "synonymous_variant")
`%notin%` <- Negate(`%in%`)
gwascatalog_AD <- gwascatalog_AD[gwascatalog_AD$CONTEXT %notin% coding, ]
which_ad <- stringr::str_detect(gwascatalog_AD$DISEASE.TRAIT, "Alzheimer")
gwascatalog_AD <- gwascatalog_AD[which_ad, ]

colnames(gwascatalog_AD)[29] <- "P"
colnames(gwascatalog_AD)[23] <- "rsid"
colnames(gwascatalog_AD)[14] <- "bp"
colnames(gwascatalog_AD)[13] <- "Chr"

gwascatalog <- gwascatalog_AD[, c("rsid", "Chr", "bp", "P")]
readr::write_csv(gwascatalog, "data/tofix_gwascatalog.csv")
gwascatalog <- readr::read_csv("data/tofix_gwascatalog.csv")

Jansen_SNP <- vroom::vroom("data/Supplementary_Table_8_Jansen.csv")
Jansen_SNP <- Jansen_SNP[Jansen_SNP$`Credible Causal`=="Yes",]
Jansen_SNP <- Jansen_SNP[Jansen_SNP$P <= 5E-8,]
colnames(Jansen_SNP)[2] <- "rsid"
Jansen <-  Jansen_SNP[, c("rsid", "Chr", "bp", "P")]

grasp_ad <- data.table::fread("data/GRASP_AD.csv")
colnames(grasp_ad)[2] <- "rsid"
colnames(grasp_ad)[8] <- "Chr"
colnames(grasp_ad)[9] <- "bp"
colnames(grasp_ad)[3] <- "P"
table(grasp_ad$FunctionalClass)
coding_grasp <- c("cds-synon", "missense")
grasp_ad <- grasp_ad[grasp_ad$FunctionalClass %notin% coding_grasp, ]
grasp <- grasp_ad[, c("rsid", "Chr", "bp", "P")]

credSNP <- rbind(gwascatalog, grasp, Jansen)
credSNP <- credSNP[!duplicated(credSNP$rsid), ]

credranges <- GRanges(credSNP$Chr, IRanges(credSNP$bp, credSNP$bp), rsid = credSNP$rsid, P=credSNP$P)
save(credranges, file ="data/AD_credibleSNP.rda")

# "http://genome.ucsc.edu/cgi-bin/hgTables"
