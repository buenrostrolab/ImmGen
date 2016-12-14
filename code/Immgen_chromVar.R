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
# Performed by Caleb Lareau, 19 November

#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
csv <-  "/Users/lareauc/Dropbox (Partners HealthCare)/Immgen_chromVAR/ImmGen_ATAC_Counts.csv"
counts <- Matrix(data.matrix(read_csv(csv)))
peakdf <- read.table("/Users/lareauc/Dropbox (Partners HealthCare)/Immgen_chromVAR/ImmGen.ATAC.master.peaks.bed")
names(peakdf) <- c("chr", "start", "end", "name")
peaks <- GRanges(peakdf)
# 
peaks <- get_gc(peaks, genome = BSgenome.Mmusculus.UCSC.mm10)

#' ## Get motifs; compute deviations. 
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
bg <- get_background_peaks(counts_mat = counts, bias = peaks)
 
# chromVAR motifs
#motifs <- get_motifs(species = "Mus musculus")
 # load our own
load("../data/cisbp_mm9_unique.08_Jun_2016.RData")
# 
motif_ix <- match_pwms(pwms, peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
motif_score <- match_pwms(pwms, peaks, genome = BSgenome.Mmusculus.UCSC.mm10, out = "scores")
deviations <- compute_deviations(counts_mat = counts, background_peaks = bg, peak_indices = motif_ix)
saveRDS(deviations, file = "../output/deviations.rds")

deviations <- readRDS("../output/deviations.rds")
variability <- compute_variability(deviations$z)
labels <- TFBSTools::name(pwms[rownames(variability)])

namesdf <- read.table("../data/libnames.txt", header = FALSE)
namesvec <- as.character(namesdf[,2])
names(namesvec) <- namesdf[,1]

#Kmer stuff
kmer_idx <- get_kmer_indices(peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
write.table(data.frame(data.matrix(kmer_idx)), file = "../data/161213_peaksXkmer.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## View plots of variable TFs and deviation heatmap
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
#plot_variability(variability, labels = labels)
boo <- which(variability$p_value_adj<10^-323)
length(boo)
write.table(data.frame(data.matrix(motif_score)), file = "../data/161213_peaksXscores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(data.matrix(motif_ix)), file = "../data/161213_peaksXmatch.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(data.matrix((deviations$z))), file = "../data/161213_motifsXsamples_zscores.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

#' ## View TF x Sample clusters
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
#z <- deviations$z
#colnames(z) <- names(namesvec[colnames(z)])
#df <- z[boo,]
#df2 <- df
#df2[df2 > 5] <- 5
#df2[df2 < -5] <- -5
#heatmaply(df, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"), limits = c(-5, 5))

#' ## View TF x Sample clusters with c(-5,5) as limits
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
#heatmaply(df2, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## View clusters of TF
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
#tfcor <- cor(t(df))
#heatmaply(tfcor, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## View clusters of samples
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
#samplecor <- cor(df)
#heatmaply(samplecor, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))
