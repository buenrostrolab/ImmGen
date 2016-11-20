#+ echo=TRUE, message=FALSE, warning=FALSE
library(readr)
library(chromVAR)
library(Matrix)
library(GenomicRanges)
source("cl.R")

if (basename(getwd()) != "code") setwd("code")

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = TRUE))

bed <- "../data/ImmGen.ATAC.master.peaks.bed"
bamfiles <- list.files("../../ImmGen_BAM", full.names = TRUE)[c(TRUE,FALSE)]
peaks <- get_peaks(bed, sort_peaks = FALSE)
counts <- get_counts_from_bams_cl(bamfiles, peaks, paired = TRUE)
colnamesLib <- sapply(strsplit(basename(bamfiles), split = "[.]"), function(s) s[1])
colnames(counts) <- colnamesLib
write.table(counts, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE, file = "../data/ImmGen_ATAC_Counts.csv")


