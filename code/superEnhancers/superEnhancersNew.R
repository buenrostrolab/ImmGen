
# Import Hide's Data
library(data.table)
library(GenomicRanges)

importSE <- function(file) {
  a <- data.frame(fread(file))
  d <- data.frame(seqnames = a$V2, start = a$V3, end = a$V4)
  g <- makeGRangesFromDataFrame(d)
  return(g)
}
Bcell <- importSE("SE.by.H3K27ac/GEO_Bcellucsc.SE")
MicroGlia <- importSE("SE.by.H3K27ac/GEO_merged.microglia.SE")
pMAC <-  importSE("SE.by.H3K27ac/GEO_merged.pMAC.SE")
spBcell <- importSE("SE.by.H3K27ac/GEO_merged.Sp.Bcell.SE")
NK <- importSE("SE.by.H3K27ac/GEO_Sp_NK.SE")
spCD4 <- importSE("SE.by.H3K27ac/Sakaguchi_CD4.SP.Tcell.Thy.SE")
CD4 <- importSE("SE.by.H3K27ac/Sakaguchi_CD4.Tcell.SE")
CD8 <- importSE("SE.by.H3K27ac/Sakaguchi_CD8.Tcell.SE")
DN <- importSE("SE.by.H3K27ac/Sakaguchi_DN.Tcell.SE")
DP <- importSE("SE.by.H3K27ac/Sakaguchi_DP.Tcell.SE")
Treg <- importSE("SE.by.H3K27ac/Sakaguchi_Treg.SE")

all <- list(Bcell, MicroGlia, pMAC, spBcell, NK, spCD4, CD4, CD8, DN, DP, Treg)
allGRL <- GRangesList(Bcell, MicroGlia, pMAC, spBcell, NK, spCD4, CD4, CD8, DN, DP, Treg)
totalRegions <- reduce(c(Bcell, MicroGlia, pMAC, spBcell, NK, spCD4, CD4, CD8, DN, DP, Treg))
majorityWins <-  GRanges(slice(coverage(allGRL), lower = 6, upper = 11, rangesOnly = TRUE))

inAll <- Reduce(subsetByOverlaps, all, maxgap = 1000)
