options(stringsAsFactors = FALSE)
source("OmniCircos2.R")

library(data.table)
library(irlba)

# TSS Peaks
d <- read.table("../../../immgen_dat/ImmGenATAC1219.peak.bed", header = FALSE)
anno <- as.character(read.table("../../data/peakAnnoVector.txt", header =TRUE)[,1])
dd <- d[anno == "TSS", ]
pddf <- data.frame(dd$V1, (dd$V2 + dd$V3)/2)

# Counts 
csv <- "../../../immgen_dat/brian/ImmGen_ATACseq_17Jan17.csv.gz"
counts <- data.matrix(data.frame(fread(input = paste0('zcat < ', csv))[anno == "TSS",]))
pcs <- irlba(counts, nv = 5, nu=0)$v
hc <- hclust(as.dist(cor(t(pcs), method = "pearson")))

# Y chromosome-less genome
gNOME <- "mm10"
pofile.s <- paste("UCSC.", gNOME, ".RData", sep="")
pofile   <- system.file("data", pofile.s, package="OmicCircos")
chr.po     <- get(load(pofile))[-21,]
maxx <- max(chr.po$seg.sum.end)
chr.po$angle.start <- 270 + (chr.po$seg.sum.start/maxx)*(630-270) + 0.5
chr.po$angle.end <- 270 + (chr.po$seg.sum.end/maxx)*(630-270)
chr.po$angle.start <- chr.po$angle.start 

pdf("immgen.tss.pdf")
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos2(R=400, cir=chr.po, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, lwd=4, col = "black")
for(i in 0:84){
  circos2(R=224+2*(i), cir=chr.po, W=2, mapping=setNames(data.frame(pddf,log2(as.numeric(counts[,colnames(counts)[hc$order[i+1]]])+1)),c("chr", "po", "vals")),
          type="immDensity", cluster=FALSE, col.bar=FALSE, lwd=0, B = FALSE)
}
# Exceptional case to show color bar
i <- 85
circos2(R=224+2*(i), cir=chr.po, W=2, mapping=setNames(data.frame(pddf, log2(as.numeric(counts[,colnames(counts)[hc$order[i+1]]])+1)), c("chr", "po", "vals")),
          type="immDensity", cluster=FALSE, col.bar=TRUE, lwd=0, B = FALSE)
circos2(R=180, cir=chr.po, W=40, mapping=setNames(data.frame(pddf[c(-1:-10),], zoo::rollmedian(matrixStats::rowVars(counts),11)), c("chr", "po", "vals")),
       col.v=3, type="l",   B=FALSE, col="orange", lwd=1, scale=FALSE)

dev.off()



