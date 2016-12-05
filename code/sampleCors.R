#' ---
#' title: "ImmGen ATAC Correlations"
#' author: "Caleb Lareau"
#' date: "December 9th, 2016"
#' ---

#+ echo=FALSE, message=FALSE, warning=FALSE
library(readr)
library(plotly)
library(heatmaply)
library(scTools)
library(BuenColors)
library(ggplot2)


if (basename(getwd()) != "code") setwd("code")

#counts <- read_csv("/Users/lareauc/Dropbox (Partners HealthCare)/Immgen_chromVAR/ImmGen_ATAC_Counts.csv.gz")
#namesdf <- read.table("../data/libnames.txt", header = FALSE)
#namesvec <- as.character(namesdf[,2])
#names(namesvec) <- namesdf[,1]
#colnames(counts) <- namesvec[colnames(counts)]
# t <- as.character(read.table("../data/labels.txt")[,1])
# outs <- counts[which(t == "outside"), ]
# TSS <- counts[which(t == "TSS"), ]
# genebody <- counts[which(t == "Genebody"), ]
# 
# cor.TSS <- cor(TSS, use="complete.obs", method="pearson")
# cor.gene <- cor(genebody, use="complete.obs", method="pearson")
# cor.outs <- cor(outs, use="complete.obs", method="pearson")
# saveRDS(cor.TSS, file = "../data/cor.TSS.rds")
# saveRDS(cor.gene, file = "../data/cor.gene.rds")
# saveRDS(cor.outs, file = "../data/cor.outs.rds")
#saveRDS(cor.all, file = "../data/cor.all.rds")
cor.outs <- readRDS("../data/cor.outs.rds")
cor.all <- readRDS("../data/cor.all.rds")
cor.TSS <- readRDS("../data/cor.TSS.rds")

#' # All Peaks
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
heatmaply(cor.all, limits = c(0.25,1),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colors = jdb_palette("solar_extra"), limits = c(0.25, 1)))

#' # Promoters
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
heatmaply(cor.TSS, limits = c(0.25,1),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colors = jdb_palette("solar_extra"), limits = c(0.25, 1)))

#' # Distal
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
heatmaply(cor.outs, limits = c(0.25,1),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colors = jdb_palette("solar_extra"), limits = c(0.25, 1)))
