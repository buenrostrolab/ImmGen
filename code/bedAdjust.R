library(GenomicRanges)

ob <- "data/120923_Immgen.bed"
gdf <- setNames(read.table(ob)[,c(1,2,3,4)], c("chr", "start", "stop", "peak"))
gdf$mid <- floor((gdf$start + gdf$stop)/2)        
gdf$start <- gdf$mid - 250
gdf$stop <- gdf$mid + 249
write.table(gdf[,1:4], row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, file = "data/021016_ImmGen_500bpPeaks.bed")

