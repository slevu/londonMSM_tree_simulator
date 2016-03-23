# devtools::install_github("emvolz-phylodynamics/phydynR")
library(phydynR)
# rm(list=ls())

###--- re-run the SA model with PhyDynR
system.time(
  source('model0_phydynr.R')
)


MH <- 10 # look up to 10 years in past for infector probs

#######################

## run model0 functions
#~ o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'euler')
o <- ode(y=y0, times=times_day, func=dydt, parms=list() , method = 'adams') 
tfgy <- .tfgy( o )
sampleTimes <- scan( file = 'sampleTimes' )
ss <- matrix( scan( file = 'sampleStates' ) , byrow=TRUE, ncol = m) 
colnames(ss) <- DEMES
# regularise
ss <- ss + 1e-4
ss = sampleStates <- ss / rowSums(ss)

## sim tree
if (F){
  print('sim tree')
  print(date())
  st.tree <- system.time( {
    tree <- sim.co.tree.fgy(tfgy, sampleTimes, sampleStates)
    })
  save(tree, file = 'phydynR-testSA0-tree.RData') 
  print(date())
} else{ # load the tree from file
  load('phydynR-testSA0-tree.RData' )
}

## cd4s & ehis
# from cori paper
#~ k =1: CD4>=500
#~ k = 2 : 350<=CD4<500.
#~ k = 3: 200<=CD4<350
#~ k = 4 : CD4<200
cd4s <- setNames( sapply( 1:nrow( tree$sampleStates), function(k){
  deme <- DEMES[ which.max(tree$sampleStates[k,] ) ]
  stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
  stage <- as.numeric( tail( strsplit(stage, '')[[1]], 1 ) )
  if (stage==1) return(1e3)
  if (stage==2) return(750)
  if (stage==3) return(400)
  if (stage==4) return(300)
  if (stage==5) return(100)
}), tree$tip.label)

ehis <- setNames( sapply( 1:nrow( tree$sampleStates), function(k){
  deme <- DEMES[ which.max(tree$sampleStates[k,] ) ]
  stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
  stage <- as.numeric( tail( strsplit(stage, '')[[1]], 1 ) ) 
  ifelse( stage==1, TRUE, FALSE)
}), tree$tip.label)

#############################
#~ incidence and prevalence
yfin <- tfgy[[4]][[length(times_day)]]
ffin <- tfgy[[2]][[length(times_day)]]
newinf <- sum(ffin[1:120, 1:120] ) * 365
plwhiv <- sum( yfin[-length(yfin)] )

# rescale tree
sampleTimes <- days2years( tree$sampleTimes )
tree$edge.length <- tree$edge.length / 365
#~ bdt <- DatedTree( tree, sampleTimes, tree$sampleStates, tol = Inf )
bdt <- DatedTree( tree, sampleTimes, tree$sampleStates, tol = .1 )

############ 
############----------------------- llllllaaaaaaaaa -----------------------------------------------------------------#################
############ 

n<- bdt$n
sampleDemes <- setNames( sapply( 1:n, function(u) DEMES[which.max( tree$sampleStates[u,])] #~ bdt$sampleDemes <- setNames( sapply( 1:n, function(u) DEMES[which.max( bdt$sampleStates
                                 st.W <- system.time( {
                                   W <- phylo.source.attribution.hiv( bdt
                                                                      , bdt$sampleTimes # must use years
                                                                      , cd4s = cd4s[bdt$tip.label] # named numeric vector, cd4 at time of sampling
                                                                      , ehi = ehis[bdt$tip.label] # named logical vector, may be NA, TRUE if patient sampl , numberPeopleLivingWithHIV = plwhiv# scalar
                                                                      , numberNewInfectionsPerYear = newinf # scalar
                                                                      , maxHeight = MH
                                                                      , res = 1e3
                                                                      , treeErrorTol = Inf
                                   ) })