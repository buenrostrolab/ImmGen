# Install CRAN Packages if necessary
install.packages("Matrix")
install.packages("readr")

# Intall Bioconductor Packages if necessary
source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
biocLite("GenomicRanges")
biocLite("BiocParallel")

# Load Libraries and CL source code
library(readr)
library(Matrix)
library(GenomicRanges)

# Register Paralleization; change 2 to 1 if multiple cores are not available
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = TRUE))  # Update this with more cores if appropriate

bed <- "/Users/lareauc/Downloads/ImmGenATAC1219.peak.bed"   # Point to new bed file with peaks
bamdir <- "/Users/lareauc/Downloads/OneDrive-2016-12-21/"  # Point to directory with bam files

bamfiles <- list.files(bamdir, full.names = TRUE, pattern = "\\.bam$")
peaks <- get_peaks(bed, sort_peaks = FALSE)
counts <- get_counts(bamfiles, peaks, paired = TRUE)  # Takes a while to execute

write.table(counts, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE, file = "../data/ImmGen_ATAC_Counts.csv")


