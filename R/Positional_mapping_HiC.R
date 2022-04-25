# Link SNPs to their putative target genes using chromatin interactions.
# Load Hi-C dataset and generate a GRange object.
hic <- read.table("data/Promoter-anchored_chromatin_loops.bed", skip=1)
colnames(hic) <- c("chr", "TSS_start", "TSS_end", "Enhancer_start", "Enhancer_end")
hicranges <- GRanges(hic$chr, IRanges(hic$TSS_start, hic$TSS_end), enhancer=hic$Enhancer_start)
olap <- findOverlaps(hicranges, promoterranges)
hicpromoter <- hicranges[queryHits(olap)]
mcols(hicpromoter) <- cbind(mcols(hicpromoter), mcols(promoterranges[subjectHits(olap)]))
hicenhancer <- GRanges(seqnames(hicpromoter), IRanges(hicpromoter$enhancer, hicpromoter$enhancer+10000),
                       gene=hicpromoter$gene)

# Overlap credible SNPs with Hi-C GRange object.
olap <- findOverlaps(credranges, hicenhancer)
credhic <- credranges[queryHits(olap)]
mcols(credhic) <- cbind(mcols(credhic), mcols(hicenhancer[subjectHits(olap)]))
