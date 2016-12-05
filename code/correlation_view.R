#' 
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
library(heatmaply)

namesdf <- read.table("../data/libnames.txt", header = FALSE)
namesvec <- as.character(namesdf[,2])
names(namesvec) <- namesdf[,1]

mat_make <- function(file){
  a <- data.matrix(read.table(file, sep = ",", header = TRUE, row.names = 1 ))
  a[is.na(a)]<- 0
  row.names(a) <- namesvec[row.names(a)]
  colnames(a) <- namesvec[colnames(a)]
  return(a)
}

#' ## Correlation No Norm
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
heatmaply(mat_make("../data/cor/correlation.no.norm.csv"),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## Only CQN
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
heatmaply(mat_make("../data/cor/correlation.only.cqn.csv"),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## Only Quantile Only
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
heatmaply(mat_make("../data/cor/correlation.only.quantile.csv"),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## CQN then Quantile
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
heatmaply(mat_make("../data/cor/correlation.cqn.then.quantile.csv"),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))

#' ## Quantile then CQN
#+ fig.width=7, fig.height=7, message = FALSE, warning = FALSE, echo=FALSE
heatmaply(mat_make("../data/cor/correlation.quantile.then.cqn.csv"),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))


