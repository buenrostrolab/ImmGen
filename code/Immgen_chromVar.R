#+ echo=TRUE, message=FALSE, warning=FALSE
library(readr)
library(chromVAR)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

#' Analysis of ImmeGen Consortium ATAC data using chromVAR
# Performed by Caleb Lareau, 24 September

#' ## Load data
#+ cache = TRUE, message=FALSE, warning=FALSE
csv <- "../data/160923_Immgen.csv.zip"
counts <- Matrix(data.matrix(read_csv(csv)[,-1]))
peakdf <- read.table("../data/120923_Immgen.bed")
names(peakdf) <- c("chr", "start", "end", "name")
peaks <- GRanges(peakdf)

peaks <- get_gc(peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
peaks <- sort(peaks)

#' ## These are the samples that didn't pass bias filtering
#+ cache = TRUE, message=FALSE, warning=FALSE
low_bias_samples <- bias_filtering(counts, peaks)
colnames(counts)[!colnames(counts) %in% names(low_bias_samples)]
counts <- counts[, low_bias_samples]

#' ### All peaks passed the filter_peaks function
#+ cache = TRUE, message=FALSE, warning=FALSE
peaks_to_keep <- filter_peaks(counts, peaks, non_overlapping = TRUE)
#counts <- counts[peaks_to_keep,]
#peaks <- peaks[peaks_to_keep] 

#' ## Get motifs; compute deviations. 
#+ cache = TRUE, message=FALSE, warning=FALSE
bg <- get_background_peaks(counts_mat = counts, bias = peaks)
motifs <- get_motifs(species = "Mus musculus")
motif_ix <- match_pwms(motifs, peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
deviations <- compute_deviations(counts_mat = counts, background_peaks = bg, peak_indices = motif_ix)
variability <- compute_variability(deviations$z)
labels <- TFBSTools::name(motifs[rownames(variability)])


#' ## View plots of variable TFs and deviation heatmap
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
plot_variability(variability, labels = labels)

boo <- which(variability$p_value_adj<0.1)
plot_deviations(deviations$z[boo,],set_names = labels[boo])

#' ## Export table of deviation scores for cells
d <- deviations$z
rownames(d) <- labels
write.table(d, file = "../output/deviationScoresImmgen.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
