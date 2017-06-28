##---- run model1.R ----
require(phydynR)
source('code/model1.R')

MH <- 20 # max height in year for SA
MU <- 1.9e-3 # mutation rate
SL <- 1e3 # seq length
DL <- 0.05 # limit of distance to record
LN <- c(0.01, 0.5) # parms log-normal relaxed clock
PID <- Sys.getpid()
OUT <- 'data/simulations2/model1-simBaseline0/'

# counterfactuals sim'ed separately, eg: 
#~ nh_wtransm <- c( 
#~ 	nh1 = 1
#~ 	,nh2 = 1
#~ 	,nh3 = 1
#~ 	,nh4 = 1
#~ 	,nh5 = 1
#~ )

#~ o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'rk4')
tfgy <- .tfgy( o )

#############################
#~ incidence and prevalence 
yfin <- tfgy[[4]][[length(times_day)]]
ffin <- tfgy[[2]][[length(times_day)]]
newinf <- sum(ffin[1:120, 1:120] ) * 365 
plwhiv <- sum( yfin[-length(yfin)] )

#~  sample time and states
sampleTimes <- scan( file = 'sampleTimes' )
ss  <- matrix( scan( file = 'sampleStates' ) , byrow=TRUE, ncol = m)
colnames(ss) <- DEMES
# regularise
ss <- ss + 1e-4
ss = sampleStates <- ss / rowSums(ss)

# modify sample times so max corredsponds to 2013 
sampleTimes <- sampleTimes + (years2days(2013) - max(sampleTimes) )

# test by downsampling
if (F)
{
n <- length( sampleTimes )
keep <- sample.int( n, size = 2e3, replace=F)
sampleTimes <- sampleTimes[keep]
sampleStates <- sampleStates[keep, ]
}


## sim tree

print('sim tree')
print(date())
st.tree <- system.time( {
#~ 	daytree <- sim.co.tree.fgy(tfgy,  sampleTimes, sampleStates)
	daytree <- rcolgem::sim.co.tree.fgy(tfgy,  sampleTimes, sampleStates)
})
print(date())



# rescale tree 
tree <- daytree
sampleTimes <- days2years( tree$sampleTimes )
tree$edge.length <- tree$edge.length / 365
bdt <- DatedTree(  tree, sampleTimes, tree$sampleStates, tol = Inf)

## compute gen distances
source('code/fn_tree2CophDist.R')
system.time( distance <- tree2CophDist2(bdt, mu = MU, seqlength = SL, dlim = DL, parmLN = LN) ) # 146s
#str(distance)
#D[[1]] # stats for all distances


n <- bdt$n

treeSampleStates <- tree$sampleStates
cd4s <- setNames( sapply( 1:nrow(treeSampleStates), function(k){
	deme <- DEMES[ which.max(treeSampleStates[k,] )  ]
	stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
	stage <- as.numeric( tail( strsplit(stage, '')[[1]], 1 ) )
	if (stage==1) return(1e3)
	if (stage==2) return(750)
	if (stage==3) return(400)
	if (stage==4) return(300)
	if (stage==5) return(100)
}), tree$tip.label)
ehis <-  setNames( sapply( 1:nrow( treeSampleStates), function(k){
	deme <- DEMES[ which.max(treeSampleStates[k,] )  ]
	stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
	stage <- as.numeric( tail( strsplit(stage, '')[[1]], 1 ) )
	ifelse( stage==1, TRUE, FALSE) 
}), tree$tip.label)
sampleDemes <- setNames( sapply( 1:n, function(u) DEMES[which.max( treeSampleStates[u,])] ), tree$tip.label )

# system.time(
#   W <- phylo.source.attribution.hiv.msm( bdt
#                                          , bdt$sampleTimes # must use years
#                                          , cd4s = cd4s[bdt$tip.label] # named numeric vector, cd4 at time of sampling 
#                                          , ehi = ehis[bdt$tip.label] # named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
#                                          , numberPeopleLivingWithHIV  = plwhiv# scalar
#                                          , numberNewInfectionsPerYear = newinf # scalar 
#                                          , maxHeight = MH 
#                                          #, res = 1e3
#                                          #, treeErrorTol = Inf
#                                          #, minEdgeLength = 1/52
#                                          , mode = 1
#   )  
# ) # 56s

#- with new sa method
system.time(
  W <- phylo.source.attribution.hiv.msm( bdt
                                         , bdt$sampleTimes # must use years
                                         , cd4s = cd4s[bdt$tip.label] # named numeric vector, cd4 at time of sampling 
                                         , ehi = ehis[bdt$tip.label] # named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
                                         , numberPeopleLivingWithHIV  = plwhiv# scalar
                                         , numberNewInfectionsPerYear = newinf # scalar 
                                         , maxHeight = MH 
                                         #, res = 1e3
                                         #, treeErrorTol = Inf
                                         #, minEdgeLength = 1/52
                                         , mode = 2
  )  
) # 698s



# str(W)
# id <- tapply(as.numeric(W$infectorProbability), W$recip, sum)
# mean(id) # [1] 0.2364833
# id2 <- tapply(as.numeric(W2$infectorProbability), W2$recip, sum)
# mean(id2) # [1] 0.4409195

save( bdt, W, cd4s, sampleDemes, plwhiv, newinf, MH, distance, file = paste(OUT, PID, '.RData', sep=''))
