require(phydynR)
source('model1.R')

PID <- Sys.getpid() 

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

D <-  cophenetic.phylo( bdt )
cat( 'mean genetic divergence if rate = .0015\n' )
print( mean(D ) * .0015 )

n<- bdt$n

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



