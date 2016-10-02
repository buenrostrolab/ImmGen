#+ echo=TRUE, message=FALSE, warning=FALSE
library(readr)
library(chromVAR)
library(Matrix)
library(GenomicRanges)

if (basename(getwd()) != "code") setwd("code")

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = TRUE))


bed <- "../data/021016_ImmGen_500bpPeaks.bed"
bamfiles <- list.files("../../../../ImmGenBam", full.names = TRUE)
peaks <- get_peaks(bed)
counts <- get_counts_from_bams_CL(bamfiles, peaks, paired = TRUE)
colnamesLib <- sapply(strsplit(basename(bamfiles), split = "[.]"), function(s) s[1])
colnames(counts) <- colnamesLib
write.table(counts, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE, file = "../data/021016_newCounts.csv")

