peaktbl <- read.table("ImmGenATAC1219.peak.bed")[,1:3]
names(peaktbl) <- c("seqnames", "start", "end")
peaks <- GenomicRanges::makeGRangesFromDataFrame(peaktbl, keep.extra.columns = TRUE)
phyloP <- "mm10.60way.phyloP60way.bw"

xx <- sapply(1:length(peaks), function(i){
    pp <- rtracklayer::import.bw(phyloP, which = peaks[i], as = "NumericList")
    if(length(pp[[1]]) == 251){
        c(i, pp[[1]])
    } else {
        c(i, rep(0, 251))   
    }
})
saveRDS(xx, file = "phyloP.raw.rds")