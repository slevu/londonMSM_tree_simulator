### source attribution on simulated trees ###
### ------------------------------------- ###
### First assume simulated tree is dated 
### as would the real tree after LSD (later try BEAST) 
### Second input dates and tree on Rcolgem
### Third analyse out-degrees
install.packages("rcolgem", repos="http://R-Forge.R-project.org")
rm(list=ls())
####---- lib ----
library(ape)

####---- load sim ----
## Load newest Rdata
l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
load(l[length(l)])
ls()

####---- look tree ----
class(tree)
str(tree)
tree
names(tree)
is.rooted(tree)

head(tree$sampleTimes)
tail(tree$sampleTimes)

tree$edge[1:1000, 1:2]
date0 <- as.Date('1979-01-01')
hist(tree$sampleTimes + date0,  "years")
hist(tree$edge.length/365)
hist(tree$heights/365)
hist(tree$sortedSampleHeights/365)
head(tree$sortedSampleHeights, 10)

##- see W0_rerunSLV.R for re-run of LSD on real data
# STFN <- '../phylo-uk/data/ExaML_result.subUKogC_noDRM.reroot_dropOG.dates'
# dates <-( read.table(STFN, skip=1, colClasses=c('character', 'numeric') ) )
# head(dates)
# sampleTimes <- scan( file = 'sampleTimes' )
# head(sampleTimes)

require(rcolgem)
require(phytools)
require(Rcpp)

## approx incidence & prevalence  London
#1. Yin Z, Brown AE, Hughes G, Nardone A, Gill ON, Delpech VC, et al. HIV  in  the  United  Kingdom  2014  Report:  data  to  end  2013 [Internet]. London: Public Health England; 2014 [cited 2015 Oct 12]. Available from: url
Y <- 43510 / 2
#2. Birrell PJ, Gill ON, Delpech VC, Brown AE, Desai S, Chadborn TR, et al. HIV incidence in men who have sex with men in England and Wales 2001–10: a nationwide population study. The Lancet Infectious Diseases. 2013 Apr;13(4):313–8.
I <-  2500 / 365 # per day

## load tree 
t <- tree
## load dates
dates <- tree$sortedSampleHeights
head(dates)
##- sample times
sts <- tree$sampleTimes 
head(sts)

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

#~ bdt <- DatedTree(t, sts, tol = .1) #slow 
#~ clade_nodes <- bdt$edge[sapply( 1:nrow(bdt$edge), function (iedge){
#~ 	a <- bdt$edge[iedge,1]
#~ 	u <- bdt$edge[iedge,2]
#~ 	if ( bdt$heights[a] > maxTMRCA & bdt$heights[u] < maxTMRCA) return(TRUE)
#~ 	FALSE
#~ }), 2]
#~ 
#~ 
#~ tree <- t
#~ sampleTimes <- sts
#~ treeErrorTol <- .1
#~ cluster_nodes <- bdt$edge[sapply( 1:nrow(bdt$edge), function (iedge){
#~ 	a <- bdt$edge[iedge,1]
#~ 	u <- bdt$edge[iedge,2]
#~ 	if ( bdt$heights[a] > maxTMRCA & bdt$heights[u] < maxTMRCA) return(TRUE)
#~ 	FALSE
#~ }), 2]
#~ cluster_nodes <- cluster_nodes[ cluster_nodes > bdt$n ]
#~ clusters <- lapply( cluster_nodes, function(u){
#~ 	t <- extract.clade( tree, u )
#~ 	DatedTree(t, sampleTimes[t$tip.label], tol = treeErrorTol)
#~ })
#~ x <- sapply( clusters, function(cl) length(cl$tip.label) )

Ws <- o2
nw <- sum( sapply( Ws, function(w) nrow(w)^2 - nrow(w) ) )
donor <- rep(NA, nw)
recipient <- rep(NA, nw)
W <- rep(NA, nw )
k <- 1
for (w in Ws ){
  for (u in 1:nrow(w)) for (v in 1:nrow(w)){
    if (u!=v){
      donor[k] <- rownames(w)[u]
      recipient[k] <- rownames(w)[v]
      W[k] <- w[u,v] 
      k <- k + 1
    }
  }
}
W_ew <- data.frame( donor, recipient,  W)


OFN0 <- 'W0.sim.RData'
OFN1 <- 'W0.sim.csv' 
save( Ws, file = OFN0 )
write.csv( W_ew, file = OFN1, quote = F , row.names=F)

