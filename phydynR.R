# devtools::install_github("emvolz-phylodynamics/phydynR")
library(phydynR)
# rm(list=ls())
# 
source('model0_phydynr.R')

## Load newest Rdata
l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
load(l[length(l)])
ls()

names(tree)
head(tree$heights)
head(tree$sampleStates)
head(tree$sampleTimes)
## approx incidence & prevalence  London
prev <- 43510 / 2
inc  <-  2500 
##- sample times
sts <- tree$sampleTimes / 365
head(sts)
## load dates
# dates <- tree$sortedSampleHeights
# head(dates)

cd4s <- setNames( sapply( 1:nrow(tree$sampleStates), function(i){
  x <- which.max( tree$sampleStates[i, ] )
  if (x ==1) return (600)
  if (x ==2 ) return (400)
  if (x==3) return (150) 
}), names( tree$sampleTimes ) )

head(cd4s)

ehi <- setNames( sapply( 1:nrow(sampleStates), function(i){
  x <- which.max( sampleStates[i, ] )
  if (x == 1) return(TRUE)
  return(FALSE)
}), names(sampleTimes ) )

## run phylo.source.attribution.hiv

# tree
# , sampleTimes # must use years
# , cd4s # named numeric vector, cd4 at time of sampling 
# , ehi # named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
# , numberPeopleLivingWithHIV # scalar
# , numberNewInfectionsPerYear # scalar 
# , maxHeight
# , res = 1e3
# , treeErrorTol = 1e-2
# 

st <- system.time(
test <- phylo.source.attribution.hiv(tree = tree,
           sampleTimes =  sts,
           cd4s = NA,
           ehi = NA,
           numberPeopleLivingWithHIV = prev,
           numberNewInfectionsPerYear = inc,
           maxHeight = 10)
)

## run it 
#~ print(date())
#~ o <- phylo.source.attribution(t, sts, f.t=I, Y.t=Y, maxTMRCA = 365*10, res = 1000, treeErrorTol = 0.1 )
#~ print(date())



sourceCpp('phylo.source.attribution.cpp')
.phylo.source.attribution <- function (tree, sampleTimes, f.t, Y.t, maxTMRCA = NULL, res = 1000, 
                                       treeErrorTol = 0.001) 
{
  if (is.null(maxTMRCA)) 
    maxTMRCA <- Inf
  if (!is.function(f.t)) {
    ft <- f.t
    f.t <- function(t) ft
  }
  if (!is.function(Y.t)) {
    Yt <- Y.t
    Y.t <- function(t) Yt
  }
  mst <- max(sampleTimes)
  sampleHeights <- mst - sampleTimes
  n <- length(sampleTimes)
  A.t <- function(t, bdt) {
    h <- mst - t
    sum(sampleHeights <= h) - sum(bdt$heights[(n + 1):(n + 
                                                         bdt$Nnode)] <= h)
  }
  d.Lambda <- function(h, Lambda, parms, ...) {
    t <- mst - h
    f <- f.t(t)
    Y <- Y.t(t)
    A <- A.t(t, parms$bdt)
    list(c(Lambda = unname(f * (Y - A)/Y^2)))
  }
  {
  bdt <- DatedTree(tree, sampleTimes, tol = treeErrorTol) # a bit slows
  heights <- seq(0, min(bdt$maxHeight, maxTMRCA), length.out = res)
  Lambda <- ode(c(Lambda = 0), times = heights, func = d.Lambda, 
                parms = list(bdt = bdt), method = "adams")[, 2]
  dh <- heights[2] - heights[1]
  Lambda.h <- function(h) Lambda[min(length(Lambda), 1 + 
                                       floor(h/dh))]
  
  ## get clusters
  cluster_nodes <- bdt$edge[sapply( 1:nrow(bdt$edge), function (iedge){
    a <- bdt$edge[iedge,1]
    u <- bdt$edge[iedge,2]
    if ( bdt$heights[a] > maxTMRCA & bdt$heights[u] < maxTMRCA) return(TRUE)
    FALSE
  }), 2]
  cluster_nodes <- cluster_nodes[ cluster_nodes > bdt$n ]
  clusters <- lapply( cluster_nodes, function(u){
    t <- extract.clade( tree, u )
    reorder.phylo(
      DatedTree(t, sampleTimes[t$tip.label], tol = treeErrorTol)
      , order = "postorder")
  })
  
  Ws <- lapply( clusters, function(cl){
    m <-matrix(0, nrow = cl$n, ncol = cl$n )
    rownames(m) = colnames(m) <- cl$tip.label
    m
  } )
  
  for (icl in 1:length(clusters))
  {
    bdt <- clusters[[icl]]
    ancs_postorder <- c()
    for (a in bdt$edge[, 1]) if (!(a %in% ancs_postorder)) 
      ancs_postorder <- c(ancs_postorder, a)
    node_daughters <- t(sapply(1:(bdt$Nnode + n), function(u) {
      i <- which(bdt$edge[, 1] == u)
      if (length(i) == 2) {
        return(bdt$edge[i, 2])
      }
      else {
        return(c(NA, NA))
      }
    }))
    PSI <- matrix(0, nrow = bdt$n + bdt$Nnode, ncol = bdt$n)
    PSI[cbind(1:bdt$n, 1:bdt$n)] <- 1
    for (inode in 1:length(ancs_postorder)) {
      a <- ancs_postorder[inode]
      {
      u <- node_daughters[a, 1]
      v <- node_daughters[a, 2]
      dpsi_au <- exp(-(Lambda.h(bdt$heights[a]) - Lambda.h(bdt$heights[u])))
      dpsi_av <- exp(-(Lambda.h(bdt$heights[a]) - Lambda.h(bdt$heights[v])))
      PSI[a, ] <- dpsi_au * PSI[u, ] + dpsi_av * PSI[v, 
                                                     ]
      utips <- which(PSI[u, ] > 0)
      vtips <- which(PSI[v, ] > 0)
      psi_a <- PSI[a, ]
      Ws[[icl]] <- updateWCpp2(Ws[[icl]], psi_a, utips, vtips)
      PSI[a, ] <- PSI[a, ]/2
      }
    }
    print( date())
    print( bdt$n )
  }
  }
  
  Ws 
}

syst <- system.time( { 
  o2 <- .phylo.source.attribution(t, sts, f.t=I, Y.t=Y, maxTMRCA = 365*10, res = 1000, treeErrorTol = 0.1 ) 
} )