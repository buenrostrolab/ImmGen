#' ---
#' title: "ImmGen ATAC Correlations"
#' author: "Caleb Lareau"
#' date: "December 9th, 2016"
#' ---

#+ echo=FALSE, message=FALSE, warning=FALSE
library(readr)
library(plotly)
library(scTools)
library(BuenColors)
library(ggplot2)
library(reshape2)

if (basename(getwd()) != "code") setwd("code")

cor.outs <- readRDS("../data/cor.outs.rds")
cor.TSS <- readRDS("../data/cor.TSS.rds")
mco <- melt(cor.outs)
mct <- melt(cor.TSS)

plotdf <- data.frame(mco, mct)
plotdf <- plotdf[plotdf$value < 1, c(1,2,3,6)]
names(plotdf) <- c("Name1", "Name2", "Distal", "Promoter")
plotdf$anno <- apply(sapply(colnames(plotdf), function(name){paste0(name, ": ", plotdf[,name])}), 1, paste, collapse = "<br>")
plotdf$Difference <- plotdf$Promoter - plotdf$Distal

#' # 
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
names(plotdf) <- c("Name1", "Name2", "Distal", "Promoter", "anno", "Difference")
plot_ly(plotdf, x = ~Promoter, y = ~Distal, mode = "markers", color = ~Difference, text = ~anno,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) 
