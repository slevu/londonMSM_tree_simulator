#########################################################
###--- use SA method from phydynR on UK bootstrap trees
###--- dated with LSD
#########################################################

# rm(list=ls())
library(ape)
library(phydynR)

 ##- list of lsd trees
 ## filename pattern from LSD changes with LSD version !
 # list.lsd.trees <- list.files(path = "data/LSD", pattern = "result.date", full.names = TRUE) # macbook
   list.lsd.trees <- list.files(path = "data/LSD", pattern = "result_newick_date", full.names = TRUE)
 head(list.lsd.trees)
 
 ##- read first LSD tree to name 
 ## sampling times and CD4s with tip.labels
 t <- read.tree( list.lsd.trees[2] )
 str(t)
 STFN <- "data/LSD/t000.dates"

 ##- parameter phydynR
 ##- Maximum height
 MH <- 10 # look up to 10 years in past for infector probs
 ##- CD4 values
 load("../phylo-uk/data/sub.RData")
 rm(s)
   ## selection of df covariates
   names(df)
   cd4s <- setNames(df$cd4, df$seqindex)[t$tip.label]
   head(cd4s)
 ##- incidence, prevalence: central scenario # todo: range of values
 ## Yin et al. 2014: 2,820 (95% CrI 1,660-4,780)
 newinf <- 2500 # c(1660, 4780)
 plwhiv <- 43150 / 2 # c(43510 / 2, 43510 / 1.5) 
 ##- sampling times
 dates <-( read.table(STFN, skip=1, 
                      colClasses=c('character', 'numeric') ) )
 head(dates)
 ##- named vector
 sampleTimes <- setNames( dates[,2], dates[,1] )[t$tip.label] 
 head(sampleTimes)
 
###--- function to apply SA over lsd tree
sa <- function(lsd_tree){
  W <- phylo.source.attribution.hiv( lsd_tree, 
          sampleTimes, # must use years
          cd4s = cd4s, # named numeric vector, cd4 at time of sampling
          ehi = NA, # named logical vector, may be NA, TRUE if patient sampl
          numberPeopleLivingWithHIV = plwhiv, # scalar
          numberNewInfectionsPerYear = newinf, # scalar
          maxHeight = MH,
          res = 1e3,
          treeErrorTol = Inf)
  return(W)
}

## loop
for (i in 1:length(list.lsd.trees)){
 tree <- read.tree(file = list.lsd.trees[i])
 W <- sa(lsd_tree = tree)
 saveRDS(W, file = paste("data/phydynR/W0_uk_", i, ".rds", sep = ''))
}

####---- reprise ----####
## First infectorprob file
list.W0 <- list.files("data/phydynR", full.names = TRUE)
## sort
list.W0 <- list.W0[order(nchar(list.W0), list.W0)]

## pick first one
 W <- readRDS(list.W0[1])
 
###--- analyses ---###
##- calculating outdegree
  ## NA values
  sum(is.na(W$infectorProbability))
out0 <- aggregate(x = list(outdegree = W$infectorProbability),
                  by = list(patient = W$donor), 
                  FUN = function(x) sum(x, na.rm = T))
head(out0)
summary(out0$outdegree)
hist(out0$outdegree)

##- add individual explanatory variates
##- selection of df covariates
y <- df[,c("seqindex","patientindex", "dob_y",
            "agediag", "cd4", "ydiag", "CHICflag")]
out <- merge(out0, y,
             by.x = "patient",
             by.y = "seqindex",
             all.x = T, sort = FALSE)

## categorize continuous variables
str(out)
age2quantile <- function(age){
if (is.na(age)) return (NA)
if (age < 27) return(1)
if (age < 33) return(2)
if (age < 40) return(3)
return(4)
}

cd4toStage <- function(cd4){
  if (is.na(cd4)) return(NA)
  if (cd4 > 700 ) return(1) #based on .9 quantile
  if (cd4 > 500 ) return(2)
  if (cd4 > 350 ) return(3)
  if (cd4 > 200 ) return(4)
  return(5)
}
out$agecl <- sapply( out[ , "agediag"] , age2quantile )
out$cd4cl <- sapply( out[ , "cd4"] , cd4toStage )
head(out)


?boxplot

boxplot(out$outdegree ~ out$agecl)
boxplot(out$outdegree ~ out$cd4cl)
boxplot(out$outdegree ~ out$CHICflag)

plot(sqrt(out$cd4), out$outdegree)
plot(out$CHICflag, out$outdegree)

############################################
###########-------- la ------############### 
############################################    

        
source("functions.R")
reg.sum.bs

model0 <- "scale(outdegree) ~ scale(agediag)"
model1 <- "scale(outdegree) ~ scale(sqrt(cd4))"
model2 <- "scale(outdegree) ~ factor(ethn.bin)"
model3 <- "scale(outdegree) ~ factor(CHICflag)"
model4 <- "scale(outdegree) ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn.bin) + factor(CHICflag)"
models <- as.list(paste0("model", 0:4))

## example
 reg.sum.bs(ls = list(out), reg = lm, model = model0) 
 #### does not work for now, see function !!!!!!!!

 ## categorize continuous variable as sims
 ## 

####---- all lm ----
lapply(models, function(x) {reg.sum.bs(ls = listUKclus, reg = lm, model = x)
})

####---- logistic ----
reg.sum.bs(ls = listUKclus, reg = glm, model = logit_model_uk, family = binomial(link = "logit"))


## from simulated trees
# ratesEqualFNS <- list.files('RData', full.names=T, path = 'simulations')
# load(ratesEqualFNS[1])
# head(W$donor)
# head(W$recip)
# head(W$infectorProbability)
           
                                 
##### from Rcolgem SA
## approx incidence & prevalence  London
#1. Yin Z, Brown AE, Hughes G, Nardone A, Gill ON, Delpech VC, et al. HIV  in  the  United  Kingdom  2014  Report:  data  to  end  2013 [Internet]. London: Public Health England; 2014 [cited 2015 Oct 12]. Available from: url
Y <- 43510 / 2
#2. Birrell PJ, Gill ON, Delpech VC, Brown AE, Desai S, Chadborn TR, et al. HIV incidence in men who have sex with men in England and Wales 2001–10: a nationwide population study. The Lancet Infectious Diseases. 2013 Apr;13(4):313–8.
I <-  2500 / 365 # per day

## load tree and dates
t <- read.tree( TFN )
#~ dates <- data.matrix( read.table(STFN, skip=1) )
dates <-( read.table(STFN, skip=1, colClasses=c('character', 'numeric') ) )
head(dates)
##- sample times
sts <- setNames( dates[,2], dates[,1] )[t$tip.label] 
head(sts)

## run it 
#~ print(date())
#~ o <- phylo.source.attribution(t, sts, f.t=I, Y.t=Y, maxTMRCA = 365*10, res = 1000, treeErrorTol = 0.1 )
#~ print(date())



sourceCpp( '~/Documents/phylo-uk/code/phylo.source.attribution.cpp')
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
      Ws[[icl]] <- updateWCpp(Ws[[icl]], psi_a, utips, vtips) # updateWCpp2 ?
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
  o <- .phylo.source.attribution(t, sts, f.t=I, Y.t=Y, maxTMRCA = 365*10, res = 1000, treeErrorTol = 0.1 ) 
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

Ws <- o
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


OFN0 <- 'W0.RData'
OFN1 <- 'W0.csv' 
save( Ws, file = OFN0 )
write.csv( W_ew, file = OFN1, quote = F , row.names=F)