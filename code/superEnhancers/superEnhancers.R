options(stringsAsFactors = FALSE)

library(data.table)
library(RcppRoll)
library(GenomeInfoDb)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)

# TSS Peaks
d <- read.table("../../../immgen_dat/ImmGenATAC1219.peak.bed", header = FALSE)
anno <- as.character(read.table("../../data/peakAnnoVector.txt", header =TRUE)[,1])
dd <- d[anno == "outside", ]
mid <- (dd$V2 + dd$V3)/2
chr_name <- dd$V1

# Counts 
csv <- "../../../immgen_dat/ATAC.density.population.mean.csv.gz"
counts <- rowSums(data.matrix(data.frame(fread(input = paste0('zcat < ', csv))[anno == "outside",-1])))

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
g$score <- as.numeric(as.character(g$score))
u.tab <- as(g, "UCSCData")
si <- Seqinfo(genome="mm10")
u.tab@seqinfo <- si
export(u.tab, format = "BigWig", con="superEnhancerScores.bw")

# Table of top loci
bedgraph$score <- as.numeric(as.character(bedgraph$score))
top <- head(bedgraph[order(bedgraph$score, decreasing = TRUE),], 100)
dout <- data.frame(paste0(top$chr, ":", as.character(as.numeric(top$start) - 10000), "-", as.character(as.numeric(top$end) + 10000)), top$score)
write.table(dout, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", file = "superEnhancerRegions.txt")

## Plots
pdf("superEnhancerRanking.pdf")
qplot(1:length(p), p) + geom_point(data = permdf, inherit.aes = FALSE, aes(rank, val, colour = "Permuted")) +
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1)) +labs(title = "Possible Super Enhancers in ImmGen", 
    x = "Rank Ordering", y = "Z Score") + theme_bw()
dev.off()

pdf("superEnhancerManhattan.pdf")
for(chr in chrs){
  p <- qplot(medi[[chr]], ((sm[[chr]]/rng[[chr]])-mean(numer))/stdev) + labs(title = chr, 
    x = "Median Peak Locus", y = "Z Score") + theme_bw()
  print(p)
}
dev.off()

