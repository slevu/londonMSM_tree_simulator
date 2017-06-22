tr <- rtree(n = 6)
plot(tr)
tr$edge.length <- tr$edge.length*20
cod <- cophenetic.phylo(tr)

fn1 <- function(m){
  require(reshape2)
  m[upper.tri(m, diag=TRUE)] <- NA
  el <- melt(m, na.rm = TRUE)
  colnames(el) <- c('ID1', 'ID2', 'distance')
  rownames(el) <- NULL
  el
}

fn2 <- function(m){
  m[upper.tri(m, diag=TRUE)] <- NA
  df <- as.data.frame(as.table(m))
  df <- df[!is.na(df$Freq), ]
  colnames(df) <- c('ID1', 'ID2', 'distance')
  rownames(df) <- NULL
  df
}

fn3 <- function(m){
  xy <- t(combn(colnames(m), 2))
  df <- data.frame(xy, dist=m[xy])
  colnames(df) <- c('ID1', 'ID2', 'distance')
  df
}

m <- cod
summary( m[lower.tri(m, diag = FALSE)] )
library(microbenchmark)
microbenchmark(fn1(cod),fn2(cod), fn3(cod))

fn1(cod)
fn3(cod)
