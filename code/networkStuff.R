#' ---
#' title: "ImmGen ATAC Plots"
#' author: "Caleb Lareau"
#' date: "December 9th, 2016"
#' ---

#+ echo=FALSE, message=FALSE, warning=FALSE
library(igraph)
library(readr)
library(plotly)
library(heatmaply)
library(scTools)
library(BuenColors)
library(ggplot2)
library(Rtsne)

if (basename(getwd()) != "code") setwd("code")

#+ echo=FALSE, message=FALSE, warning=FALSE
zscore <- read.table("/Users/lareauc/Dropbox (Partners HealthCare)/Immgen_chromVAR/161119_motifsXsamples_zscores.tsv.gz", row.names = 1)
zsn <- sweep(zscore, 2, colMeans(zscore), FUN="/")

namesdf <- read.table("../data/libnames.txt", header = FALSE)
namesvec <- as.character(namesdf[,2])
names(namesvec) <- namesdf[,1]
colnames(zsn) <- namesvec[colnames(zsn)]
colnames(zscore) <- namesvec[colnames(zscore)]

#sc <- spRing(t(data.matrix(zsn)), k = 4)
#saveRDS(sc, file = "../data/springCoords.rds")
sc <- readRDS("../data/springCoords.rds")
#tc <- Rtsne(t(data.matrix(zsn)), pca = FALSE)$Y

classes <- sapply(strsplit(colnames(zsn), split = "[.]"), function(s) s[[1]])

df <- data.frame(class = classes, sampleName = colnames(zsn))
anno <- apply(sapply(colnames(df), function(name){paste0(name, ": ", df[,name])}), 1, paste, collapse = "<br>")
plot.df <- data.frame(sc[,1], sc[,2], classes, Annotation = anno)

#' # Per-sample class membership
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
colnames(plot.df) <- c("x", "y","Color" , "Annotation")
plot_ly(plot.df, x = ~x, y = ~y, mode = "markers", color = ~Color, text = ~Annotation,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) %>%
  add_annotations(x = c(-135, 120, -125, 55, -65, 60 , 70, -50, 30, -60),
                  y = c(50, 25, -50, 55,     140, -60,160, -50, -30, -140),
                  text = c("Bcells", "Bcells", "NK", "Pre-T", "TGD", "Dendritic", "ILC3", "BECs", "Tcells", "Progenitor"),
                  xref = "x",
                  yref = "y",
                  showarrow = FALSE,
                  ax = 20,
                  ay = -40,
                  font = list(size = 16))


#' #Gata1
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
col <- as.numeric(zscore[which(grepl("Gata1", row.names(zscore))),])
df <- data.frame(class = classes, sampleName = colnames(zscore), Gata1 = col)
anno <- apply(sapply(colnames(df), function(name){paste0(name, ": ", df[,name])}), 1, paste, collapse = "<br>")
plot.df <- data.frame(sc[,1], sc[,2], col, Annotation = anno)

colnames(plot.df) <- c("x", "y", "Gata1" , "Annotation")
plot_ly(plot.df, x = ~x, y = ~y, mode = "markers", color = ~Gata1, text = ~Annotation,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) %>%
  add_annotations(x = c(-135, 120, -125, 55, -65, 60 , 70, -50, 30, -60),
                  y = c(50, 25, -50, 55,     140, -60,160, -50, -30, -140),
                  text = c("Bcells", "Bcells", "NK", "Pre-T", "TGD", "Dendritic", "ILC3", "BECs", "Tcells", "Progenitor"),
                  xref = "x",
                  yref = "y",
                  showarrow = FALSE,
                  ax = 20,
                  ay = -40,
                  font = list(size = 16))

#' # Rora
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
col <- as.numeric(zscore[which(grepl("Rora", row.names(zscore))),])
df <- data.frame(class = classes, sampleName = colnames(zscore), Rora = col)
anno <- apply(sapply(colnames(df), function(name){paste0(name, ": ", df[,name])}), 1, paste, collapse = "<br>")
plot.df <- data.frame(sc[,1], sc[,2], col, Annotation = anno)

colnames(plot.df) <- c("x", "y", "Rora" , "Annotation")
plot_ly(plot.df, x = ~x, y = ~y, mode = "markers", color = ~Rora, text = ~Annotation,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) %>%
  add_annotations(x = c(-135, 120, -125, 55, -65, 60 , 70, -50, 30, -60),
                  y = c(50, 25, -50, 55,     140, -60,160, -50, -30, -140),
                  text = c("Bcells", "Bcells", "NK", "Pre-T", "TGD", "Dendritic", "ILC3", "BECs", "Tcells", "Progenitor"),
                  xref = "x",
                  yref = "y",
                  showarrow = FALSE,
                  ax = 20,
                  ay = -40,
                  font = list(size = 16))
#' # Tbx1
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
col <- as.numeric(zscore[which(grepl("Tbx1", row.names(zscore)))[1],])
df <- data.frame(class = classes, sampleName = colnames(zscore), Tbx1 = col)
anno <- apply(sapply(colnames(df), function(name){paste0(name, ": ", df[,name])}), 1, paste, collapse = "<br>")
plot.df <- data.frame(sc[,1], sc[,2], col, Annotation = anno)

colnames(plot.df) <- c("x", "y", "Tbx1" , "Annotation")
plot_ly(plot.df, x = ~x, y = ~y, mode = "markers", color = ~Tbx1, text = ~Annotation,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) %>%
  add_annotations(x = c(-135, 120, -125, 55, -65, 60 , 70, -50, 30, -60),
                  y = c(50, 25, -50, 55,     140, -60,160, -50, -30, -140),
                  text = c("Bcells", "Bcells", "NK", "Pre-T", "TGD", "Dendritic", "ILC3", "BECs", "Tcells", "Progenitor"),
                  xref = "x",
                  yref = "y",
                  showarrow = FALSE,
                  ax = 20,
                  ay = -40,
                  font = list(size = 16))

#' # Pax5
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
col <- as.numeric(zscore[which(grepl("Pax5", row.names(zscore)))[1],])
df <- data.frame(class = classes, sampleName = colnames(zscore), Pax5 = col)
anno <- apply(sapply(colnames(df), function(name){paste0(name, ": ", df[,name])}), 1, paste, collapse = "<br>")
plot.df <- data.frame(sc[,1], sc[,2], col, Annotation = anno)

colnames(plot.df) <- c("x", "y", "Pax5" , "Annotation")
plot_ly(plot.df, x = ~x, y = ~y, mode = "markers", color = ~Pax5, text = ~Annotation,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) %>%
  add_annotations(x = c(-135, 120, -125, 55, -65, 60 , 70, -50, 30, -60),
                  y = c(50, 25, -50, 55,     140, -60,160, -50, -30, -140),
                  text = c("Bcells", "Bcells", "NK", "Pre-T", "TGD", "Dendritic", "ILC3", "BECs", "Tcells", "Progenitor"),
                  xref = "x",
                  yref = "y",
                  showarrow = FALSE,
                  ax = 20,
                  ay = -40,
                  font = list(size = 16))
