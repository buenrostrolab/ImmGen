library(reshape2)
library(ggplot2)
library(BuenColors)
library(SummarizedExperiment)

p05dev <- readRDS("../../output/p05dev.rds")
p05devqNORM <- readRDS("../../output/p05dev_qNORM.rds")
p05devNORM <- readRDS("../../output/p05devNORM.rds")

none <- melt(assays(p05dev)$z)
qnorm <- melt(assays(p05devqNORM)$z)
cqnnorm <- melt(assays(p05devNORM)$z)

df <- data.frame(None = none$value, Quant = qnorm$value, CQN = cqnnorm$value)[sample(1:dim(none)[1], 12222),]
df2 <- df[rowSums(is.na(df)) < 1,]
density <- sapply(1:dim(df2)[1], function(i){
  sum(((df2[i,1] - df2[,1])^2 + (df2[i,2] - df2[,2])^2) < 0.1)
})
df2$logDensity <- log10(density)

ggplot(df2, aes( x = Quant,  y = CQN, color = logDensity)) + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1)) +labs(title = "Sample-Motif TF Score Threshold Comparison", 
    x = "log2 TF Score, Quant", y = "log2 TF Score, CQN", colour = "Log Density") + 
    geom_point(size = 0.01) + theme_bw() + scale_color_gradientn(colors = jdb_palette("samba_night"))