options(stringsAsFactors = FALSE)

library(data.table)
library(RcppRoll)
library(GenomeInfoDb)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)

# TSS Peaks
d <- read.table("../../../immgen_dat/ImmGenATAC1219.peak.bed", header = FALSE)
anno <- as.character(read.table("../../data/peakAnnoVector.txt", header =TRUE)[,1])
#dd <- d[anno == "outside", ]
dd <- d
mid <- (dd$V2 + dd$V3)/2
chr_name <- dd$V1
gene_anno <- as.character(read.table("geneNames.txt", header =TRUE)[,1])
#gene_anno <- gene_anno[anno=="outside"]

gdf <- makeGRangesFromDataFrame(data.frame(chr = dd$V1, start = dd$V2, end = dd$V3, gene = gene_anno), keep.extra.columns = TRUE)

# Counts 
csv <- "../../../immgen_dat/ATAC.density.population.mean.csv.gz"
#counts <- rowSums(data.matrix(data.frame(fread(input = paste0('zcat < ', csv))[anno == "outside",-1])))
counts <- rowSums(data.matrix(data.frame(fread(input = paste0('zcat < ', csv))[,-1])))

chrs <- paste0("chr", c(as.character(1:19), "X"))
nn <- 7

# Rolling sum of counts
sm <- sapply(chrs, function(chr){
 roll_sum(counts[chr == chr_name], n = nn)
})

pc <- sapply(chrs, function(chr){
  roll_sum(sample(counts, replace = TRUE)[chr == chr_name], n = nn)
})

rng <- sapply(chrs, function(chr){
 roll_max(mid[chr == chr_name], n = nn) - roll_min(mid[chr == chr_name], n = nn)
})

medi <- sapply(chrs, function(chr){
 roll_max(mid[chr == chr_name], n = nn)
})

stdev <- sqrt(var(unname(unlist(pc))/unname(unlist(rng))))
numer <- unname(unlist(sm))/unname(unlist(rng))
p <- sort(numer- mean(numer))/stdev
permdf <- data.frame(1:length(p), sort(unname(unlist(pc))/unname(unlist(rng)) - mean(numer))/stdev)
names(permdf) <- c("rank", "val")

# Export to bedgraph
xx <- lapply(chrs, function(chr){
  data.frame(chr = chr, start = as.character(medi[[chr]] - 124), end = as.character(medi[[chr]] + 124), score = as.character(((sm[[chr]]/rng[[chr]])-mean(numer))/stdev))
})
bedgraph <- data.frame(rbindlist(xx, use.names=FALSE, fill=FALSE, idcol=NULL))
g <- GRanges(bedgraph)

# Export bigwig
# g$score <- as.numeric(as.character(g$score))
# u.tab <- as(g, "UCSCData")
# si <- Seqinfo(genome="mm10")
# u.tab@seqinfo <- si
# export(u.tab, format = "BigWig", con="superEnhancerScores.bw")

# Table of top loci
bedgraph$score <- as.numeric(as.character(bedgraph$score))
bedgraph$gene <- gene_anno[data.frame(findOverlaps(g, gdf))$subjectHits]
top <- head(bedgraph[order(bedgraph$score, decreasing = TRUE),], 100)
dout <- data.frame(paste0(top$chr, ":", as.character(as.numeric(top$start) - 10000), "-", as.character(as.numeric(top$end) + 10000)), top$score)
write.table(dout, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", file = "superEnhancerRegions.txt")

# Make fancier DF with names
bg <- bedgraph[order(bedgraph$score, decreasing = TRUE),]
bg <- data.frame(rank = rev(1:dim(bg)[1]), bg)

# Get Super Enhancers we like from Previous
ov <- findOverlaps(makeGRangesFromDataFrame(bg), totalRegions)
bg$Overlap <- 1:dim(bg)[1] %in% queryHits(ov)
summary(bg[bg$Overlap,"score"])
sum(head(bg$Overlap, 100))
length(which(1:dim(bg)[1] %in% queryHits(ov)))

library(ROCR)
pp <- prediction(bg$score, bg$Overlap)
ppp <- performance(pp, 'tpr','fpr')
jpeg("se_ROC.jpg")
plot(ppp,main="Super Enhancer Prediction - ROC Curves",col="blue")
abline(0,1,col="grey")
dev.off()
performance(pp, "auc")

## Plots
pdf("superEnhancerOverlap.pdf")
bg$gene[20:dim(bg)[1]] <- ""
ggplot(data = bg, aes(rank, score, color = Overlap)) + geom_point() + scale_colour_manual(values=c("black", "firebrick")) +
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1)) +labs(title = "Possible Super Enhancers in ImmGen", 
    x = "Rank Ordering", y = "Z Score") + theme_bw()
dev.off()

pdf("superEnhancerRanking.pdf")
bg$gene[20:dim(bg)[1]] <- ""
ggplot(data = bg, aes(rank, score)) + geom_point() + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1)) +labs(title = "Possible Super Enhancers in ImmGen", 
    x = "Rank Ordering", y = "Z Score") + theme_bw() + 
    geom_text_repel(data = bg, aes(rank, score, label = gene), inherit.aes = FALSE)
#+  geom_point(data = permdf, inherit.aes = FALSE, aes(rank, val, colour = "Permuted"))
dev.off()



pdf("superEnhancerManhattan.pdf")
for(chr in chrs){
  p <- qplot(medi[[chr]], ((sm[[chr]]/rng[[chr]])-mean(numer))/stdev) + labs(title = chr, 
    x = "Median Peak Locus", y = "Z Score") + theme_bw()
  print(p)
}
dev.off()

