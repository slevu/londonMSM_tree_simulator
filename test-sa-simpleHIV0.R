#~ In this script:
# simulate a tree using simple hiv model
# compute inf probs using true model
# compute inf probs using approx hiv model, using cd4s & ehi indicator 
require(phydynR)

INFECTEDNAMES <- c('I0', 'I1', 'I2')

births <- rbind( 
  c('beta0 * S * I0 / (S + I0 + I1 + I2)', '0', '0'), 
  c('beta1 * S * I1 / (S + I0 + I1 + I2)', '0', '0'), 
  c('beta2 * S * I2 / (S + I0 + I1 + I2)', '0', '0')
)
rownames(births)=colnames(births)<- INFECTEDNAMES

migrations <- rbind(
  c('0', 'gamma0 * I0', '0'),
  c('0', '0', 'gamma1 * I1'),
  c('0', '0', '0')
)
rownames(migrations)=colnames(migrations) <- INFECTEDNAMES

deaths <- c(
 'mu*I0'
 , 'mu*I1'
 , 'mu*I2 + gamma2 * I2'
)
names(deaths) <- INFECTEDNAMES

nonDemeDynamics <- c( S = '-mu*S + mu*(S + I0 + I1 + I2) -
   S * (beta0*I0+beta1*I1+beta2*I2) / (S + I0 + I1 + I2)'
)



demo.model <- build.demographic.process( 
  births
  , nonDemeDynamics
  , migrations=migrations
  , deaths=deaths
  , parameterNames = c(
     'beta0'
     , 'beta1'
     , 'beta2'
     , 'gamma0'
     , 'gamma1'
     , 'gamma2'
     , 'mu')
  , rcpp = TRUE
  , sde=FALSE
)


theta <- c(gamma0 = 1
 , gamma1 = 1/7
 , gamma2 = 1/2
 , mu = 1/30
 , beta0 = 12./10
 , beta1=3./100
 , beta2=9./100
)
t0 <- 0
t1 <- 50
x0 <- c(S = 999, I0 = 1, I1 =.1, I2 = .1)


n <- 150 
sampleTimes <- setNames( seq( 20, 25, length.out = n),  paste(sep='.', 't', 1:n) )
sampleStates <- t(rmultinom( n, size = 1, prob = c(.025, .9, .075) ))

tree <- sim.co.tree (theta, demo.model, x0, t0, sampleTimes, sampleStates, res = 1e3)

if (T)
{
W <- phylo.source.attribution.multiDeme.model( tree
  , sampleTimes
  , sampleStates
  , maxHeight = 20
  , theta, demo.model, x0, t0
  , res = 1e3
  , treeErrorTol = 1e-2
  , integrationMethod='adams'
) 
}


## This version is designed for general purpose HIV epidemic using cd4 & EHI indicator
if (F)
{
mh <- 10 # only compute infector probs for TMRCAs up to this many years in the past
cd4s <- setNames( sapply( 1:nrow(sampleStates), function(i){
	x <- which.max( sampleStates[i, ] )
	if (x ==1) return (600)
	if (x ==2 ) return (400)
	if (x==3) return (150) 
}), names( sampleTimes ) )
ehi <- setNames( sapply( 1:nrow(sampleStates), function(i){
	x <- which.max( sampleStates[i, ] )
	if (x == 1) return(TRUE)
	return(FALSE)
}), names(sampleTimes ) )
tfgy <- demo.model(theta, x0, t0, max(sampleTimes), res = 1e3, integrationMethod='adams')
i <- which(tfgy[[1]] > (max(sampleTimes) - mh ) )
plwhiv <- sapply(i , function(ii) sum(tfgy[[4]][[ii]]) )
# harm mean
plwhiv <- 1 / (mean( 1/ plwhiv))
newinf <- sapply(i, function(ii) sum( tfgy[[2]][[ii]] ) )
#  mean
newinf <- mean( newinf)# prod(newinf)^(1/length(newinf))
{
	W <- phylo.source.attribution.hiv( tree # dated tree
	  , sampleTimes # *must* use years
	  , cd4s = cd4s # named numeric vector, cd4 at time of sampling 
	  , ehi = ehi # named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
	  , numberPeopleLivingWithHIV  = plwhiv# scalar, assumed constant
	  , numberNewInfectionsPerYear = newinf # scalar, assumed constant
	  , maxHeight = mh
	  , res = 1e3
	  , treeErrorTol = Inf
	) 
}
}


## plot tree distance vs infector prob 
if (T)
{
D <- cophenetic( tree )
tdists <- sapply( 1:length( W$donor ), function(i) {
	max( D[W$donor[i], W$recip[i]], D[W$donor[i], W$recip[i]] )
})
#~ plot( tdists, W$infectorProbability , log = 'xy')
plot( tdists, W$infectorProbability , log = 'y', ylim = c( 1e-5, 1))
}

#~ show.demographic.process( demo.model, theta, x0, t0,50, res = 1000, integrationMethod = "adams")
#~ tfgy <- demo.model(theta, x0, t0, max(sampleTimes), res = 1e3, integrationMethod='adams')
indegree <- sapply( tree$tip.label, function(tl) sum(W$infectorProbability[which(W$recip==tl)] ) )
outdegree <- sapply( tree$tip.label, function(tl) sum(W$infectorProbability[which(W$donor==tl)] ) )
print(mean(indegree))
print(mean(outdegree))
X11() ; plot( ecdf( indegree ))
X11(); plot(indegree, outdegree )
