#+ echo=TRUE, message=FALSE, warning=FALSE
library(readr)
library(chromVAR)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(plotly)
library(heatmaply)

if (basename(getwd()) != "code") setwd("code")

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

# chromVAR motifs
#motifs <- get_motifs(species = "Mus musculus")
# load our own
load("../data/cisbp_mm9_unique.08_Jun_2016.RData")

#motif_ix <- match_pwms(pwms, peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
#deviations <- compute_deviations(counts_mat = counts, background_peaks = bg, peak_indices = motif_ix)
#saveRDS(deviations, file = "../output/deviations.rds")
deviations <- readRDS("../output/deviations.rds")
variability <- compute_variability(deviations$z)
labels <- TFBSTools::name(pwms[rownames(variability)])


#' ## View plots of variable TFs and deviation heatmap
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
plot_variability(variability, labels = labels)

boo <- which(variability$p_value_adj<0.1)
#plot_deviations(deviations$z[boo,],set_names = labels[boo])

#' ## View clusters of TF
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
df <- deviations$z[boo,]
tfcor <- cor(t(df))
heatmaply(tfcor)

#' ## View clusters of samples
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
samplecor <- cor(df)
heatmaply(samplecor)
