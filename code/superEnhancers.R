options(stringsAsFactors = FALSE)

library(data.table)
library(RcppRoll)

# TSS Peaks
d <- read.table("../../immgen_dat/ImmGenATAC1219.peak.bed", header = FALSE)
anno <- as.character(read.table("../data/peakAnnoVector.txt", header =TRUE)[,1])
dd <- d[anno == "outside", ]
mid <- (dd$V2 + dd$V3)/2
chr_name <- dd$V1

# Counts 
csv <- "../../immgen_dat/ATAC.density.population.mean.csv.gz"
counts <- rowSums(data.matrix(data.frame(fread(input = paste0('zcat < ', csv))[anno == "outside",-1])))

chrs <- paste0("chr", c(as.character(1:19), "X"))


nn <- 15

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