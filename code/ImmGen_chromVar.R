#+ echo=FALSE, message=FALSE, warning=FALSE
library(data.table)
library(chromVAR)
library(SummarizedExperiment)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(plotly)
library(heatmaply)

if (basename(getwd()) != "code") setwd("code")

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

#' Analysis of ImmGen Consortium ATAC data using chromVAR
# Performed by Caleb Lareau, 10 January

#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
immgenRDS <- "../output/immgenRDS.rds"
bgpeaks <- "../output/backgroundPeaks.rds"
expectations <- "../output/expectations.rds"

if(!file.exists(immgenRDS) | !file.exists(bgpeaks) | !file.exists(expectations)){
    csv <- "../../immgen_dat/total.cnt.table.1226.csv.gz"
    counts <- fread(input = paste0('zcat < ', csv))
    libnames <- counts[1,]
    libnames <- libnames[,-1]
    counts <- counts[-1,][,-1]
    counts <- Matrix(data.matrix(counts))
    
    peakdf <- data.frame(fread("zcat < ../../immgen_dat/ImmGenATAC1219.peak.bed.gz"))
    names(peakdf) <- c("chr", "start", "end", "name")
    peaks <- GRanges(peakdf)
    immgen <- SummarizedExperiment(rowRanges = peaks, colData = t(libnames), assays = list(counts = counts))
    immgen <- add_gc_bias(immgen, genome = BSgenome.Mmusculus.UCSC.mm10)
    saveRDS(immgen, file = immgenRDS)
    
    bg <- get_background_peaks(immgen)
    saveRDS(bg, file = bgpeaks)
    
    expectation <- compute_expectations(immgen)
    saveRDS(expectation, file = expectations)
}

immgen <- readRDS(immgenRDS)
bg <- readRDS(bgpeaks)
expectation <- readRDS(expectations)

#' ## Get motifs; compute deviations. 
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
if(FALSE){
  fimo <-  "../../immgen_dat/fimo_cisBP_jmotifs_p4.txt.gz"
  fimotab <- fread(input = paste0('zcat < ', fimo))
  peakidx <- as.integer(unname(sapply(data.frame(fimotab[,1])[,1], function(l) strsplit(l, split = "_")[[1]][2])))
  motifidx <- as.factor(data.frame(fimotab[,2])[,1])
  motiflvl <- levels(motifidx)
  p04match <- SummarizedExperiment(assays=list(matches = sparseMatrix(i = peakidx, j = as.numeric(motifidx))),
                                   colData = DataFrame(motif = motiflvl, row.names = motiflvl), rowRanges = immgen@rowRanges)
  saveRDS(p04match, file = "../output/match04se.rds")
  
  dev04 <- compute_deviations(immgen, annotations = p04match, background_peaks = bg, expectation = expectation)
  saveRDS(dev04, file="../output/p04dev.rds")
  
  idx <- as.numeric(data.matrix(fimotab[,4])) < 1*10^(-5)
  p05match <- SummarizedExperiment(assays=list(matches = sparseMatrix(i = peakidx[idx], j = as.numeric(motifidx)[idx])),
                                   colData = DataFrame(motif = motiflvl, row.names = motiflvl), rowRanges = immgen@rowRanges)
  dev05 <- compute_deviations(immgen, annotations = p05match, background_peaks = bg, expectation = expectation)
  saveRDS(dev05, file="../output/p05dev.rds")
}

# match object
p04match <- readRDS("../output/match04se.rds")

dev04 <- readRDS("../output/p04dev.rds")
dev05 <- readRDS("../output/p05dev.rds")
md04 <- melt(assays(dev04)[["z"]])
md05 <- melt(assays(dev05)[["z"]])

# Make Density plot using JDB method
mda <- data.frame(md04, md05)
mda <- mda[complete.cases(mda),]
X <- mda[,3]
Y <- mda[,6]
density <- pbapply::pbsapply(1:length(X), function(i){
  sum(((X[i] - X)^2 + (Y[i] - Y)^2) < 1)
})
cd <- data.frame(x = X, y = Y, density = density)
saveRDS(cd, "../output/TFscoreDesnity.rds")

cd <- readRDS("../output/TFscoreDesnity.rds")
ggplot(cd[sample(nrow(cd)),], aes( x = x,  y = y, color = density)) + 
    geom_point() + theme_bw() + scale_fill_viridis()



