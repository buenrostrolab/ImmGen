library(plotly)
library(igraph)

# Define cells and relative ordering on the plot
index <-c(             1,              2,               3,             4,              5)
Group <-c(          "SP",           "SP",           "SP",           "SP",           "SP")
cell <- c("LTHSC.34-.BM", "LTHSC.34+.BM", "STHSC.150-.BM", "MMP3.48+.BM", "MMP4.135+.BM")
X <-    c(             0,              0,               0,            -2,              2)
Y <-    c(            10,              8,               6,             4,              4)


edf <- data.frame(
  from = c(1, 2, 3, 3   ),
  to   = c(2, 3, 4, 5   ) 
)


df <- data.frame(X, Y, cell)
makeEdge <- function(row1, row2, df){
  return(c(df[row1,"X"], df[row1,"Y"], df[row2,"X"], df[row2,"Y"]))
}
network <- plot_ly(x = ~df$X, y = ~df$Y, text = df$cell,
                   type="scatter", mode="makers", marker = list(size = 30)) %>%
  add_text(textposition = "right")


edge_shapes <- list()
for(i in 1:dim(edf)[1]) {

  e <- makeEdge(edf[i,1], edf[i,2], df)
  print(e)
  edge_shape = list(
    type = "line",
    line = list(color = "#030303", width = 1),
    x0 = e[1],
    y0 = e[2] - 0.2,
    x1 = e[3],
    y1 = e[4] + 0.2
  )
  edge_shapes[[i]] <- edge_shape
}
axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)


p <- layout(
  network,
  title = 'ImmGen Lineage',
  shapes = edge_shapes,
  xaxis = axis,
  yaxis = axis
)

p