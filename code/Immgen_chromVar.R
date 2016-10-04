#+ echo=FALSE, message=FALSE, warning=FALSE
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

#' Analysis of ImmGen Consortium ATAC data using chromVAR
# Performed by Caleb Lareau, 2 October

#' ## Load data
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
csv <-  "../data/021016_newCounts.csv.zip"
counts <- Matrix(data.matrix(read_csv(csv)))
peakdf <- read.table("../data/021016_ImmGen_500bpPeaks.bed")
names(peakdf) <- c("chr", "start", "end", "name")
peaks <- GRanges(peakdf)

peaks <- get_gc(peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
praw <- peaks
peaks <- sort(peaks)
counts <- counts[match(peaks,praw),]

#' ## Update sample names from lib -> cell name
#+ echo = FALSE
namesdf <- read.table("../data/libnames.txt", header = TRUE)
namesvec <- as.character(namesdf[,2])
names(namesvec) <- namesdf[,1]
colnames(counts) <- unname(namesvec[colnames(counts)])

#' ## These are the samples that didn't pass bias filtering
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
low_bias_samples <- bias_filtering(counts, peaks)
colnames(counts)[!colnames(counts) %in% names(low_bias_samples)]
counts <- counts[, low_bias_samples]

#' ### All peaks passed the filter_peaks function
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
peaks_to_keep <- filter_peaks(counts, peaks, non_overlapping = TRUE)
counts <- counts[peaks_to_keep,]
peaks <- peaks[peaks_to_keep] 

#' ## Get motifs; compute deviations. 
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
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
boo <- which(variability$p_value_adj<0.0001)

#' ## View TF x Sample clusters
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
df <- deviations$z[boo,]
df2 <- df
df2[df2 > 5] <- 5
df2[df2 < -5] <- -5
heatmaply(df, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## View TF x Sample clusters with c(-5,5) as limits
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
heatmaply(df2, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## View clusters of TF
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
tfcor <- cor(t(df))
heatmaply(tfcor, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## View clusters of samples
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
samplecor <- cor(df)
heatmaply(samplecor, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))
