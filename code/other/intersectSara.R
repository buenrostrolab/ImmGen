library(data.table)
library(Matrix)
library(ggplot2)
library(reshape2)


if (basename(getwd()) != "code") setwd("code")

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

chromVar <- Matrix(data.matrix(fread("/Users/lareauc/Dropbox\ (Partners\ HealthCare)/Immgen_chromVAR/161213_peaksXmatch.tsv")))
sara <- data.frame(fread("/Users/lareauc/Dropbox (Partners HealthCare)/Immgen_chromVAR/cisBP_fimo_minP_t6.txt"))

tfbridge <- data.frame(fread("../data/TF_Information.txt"))

# Loop over all unique factors in Sara's data; find match in TF_Information file (grab first if multiple). 
# If no direct match can easily be found, return missing designation
oo <- sapply(unique(sara$V2), function(fac){
  ens <- tfbridge[tfbridge[,4]==fac,6][1]
  if(!is.na(ens)){
    peaks <- chromVar[,grepl(ens, colnames(chromVar)), drop = FALSE]
    sarahits <- as.numeric(unlist(strsplit(sara[sara[,2] == fac,1], split = "_"))[c(FALSE,TRUE)])
    chromVarhits <- summary(peaks)$i
    c(length(setdiff(sarahits, chromVarhits)), length(intersect(sarahits, chromVarhits)), fac, ens)
  } else {
    c(-9,-9, fac, "none")
  }
})

dfbig <- data.frame(t(oo))
colnames(dfbig) <- c("InFimoNotChromVar", "Intersection", "Fimo", "ENS")
dfplot <- dfbig[dfbig[,4] != "none",]
dfplot$PerfectMatch <- dfplot[,1] == 0 
dfplot$InFimoNotChromVar <- as.numeric(dfplot$InFimoNotChromVar)
dfplot$Intersection <- as.numeric(dfplot$Intersection)

ggplot(dfplot, aes(x = Intersection, y = InFimoNotChromVar, col = PerfectMatch)) + geom_point() + theme_bw()


