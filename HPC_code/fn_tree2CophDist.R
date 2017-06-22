tree2CophDist <- function(yeartree, mu = .0015, seqlength = 1e3, dlim = NULL){
  require(ape)
  require(reshape2)
  set.seed(123)
  # (year * susbt/site/year * n site) = number of substitutions 
  yeartree$edge.length <- rpois(length(yeartree$edge.length), yeartree$edge.length * mu * seqlength ) / seqlength 
  m <- cophenetic.phylo(yeartree)
  ###--- keep only the lower triangle by 
  ## filling upper with NA
  m[upper.tri(m, diag=TRUE)] <- NA
  ###--- create edge list if not there
  el <- melt(m, na.rm = TRUE)
  colnames(el) <- c('ID1', 'ID2', 'distance')
  if(!is.null(dlim)){
	  el <- el[el$distance < dlim,]
  }
  el
}

tree2CophDist2 <- function(yeartree, mu = .0015, seqlength = 1e3, dlim = NULL){
  require(ape)
  set.seed(123)
  # (year * susbt/site/year * n site) = number of substitutions 
  yeartree$edge.length <- rpois(length(yeartree$edge.length), yeartree$edge.length * mu * seqlength ) / seqlength 
  m <- cophenetic.phylo(yeartree)
  ##- stats on all dist
  stats <- summary( m[lower.tri(m, diag = FALSE)] )
  ##- list of unique pairs (unordered) n(n-1)/2
  xy <- t(combn(colnames(m), 2))
  ##- 3 columns df
  el <- data.frame(xy, dist = m[xy])
  colnames(el) <- c('ID1', 'ID2', 'distance')
  if(!is.null(dlim)){
	  el <- el[el$distance < dlim,]
  }
  return(list(stats = stats, dlim = dlim, D = el ))
}