# characterise 1) assort in baseline
#~ 2) clustering assort in baseline
#~ 3) out degree by age; try linreg adjusting for stage as well
THRESHOLD_YRS <- Inf #10 # only count donors sampled within this many years of the most recent sample (control for cohort effects)
CL_THRESHOLD <- .015

require(phydynR)
source('model0.R')
require(lattice)
lattice.options(default.theme = standard.theme(color=F))

ratesBaselineFNS <- list.files('RData', full.names=T, path = 'model0-simulateBaseline0')

deme2age <- function(deme){as.numeric( substr(regmatches( deme , regexpr( '\\.age[0-9]', deme )), 5,5) ) }
deme2stage <- function(deme){as.numeric( substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6,6) ) }

fn2ageMatrix_table <- function(fn)
{#~ fn <- ratesBaselineFNS[1]{
	load(fn)
	samplesInCohort <- names(daytree$sampleTimes[daytree$sampleTimes > (max(daytree$sampleTimes)-THRESHOLD_YRS*365) ] )
	samplesInCohort_donors <- intersect( samplesInCohort, W$donor)
	
	i <- which( (W$donor %in% samplesInCohort) & (W$recip %in% samplesInCohort ) )
	WW <- list( donor = W$donor[i] , recip = W$recip[i], infectorProbability = W$infectorProbability[i] )
	
	#system.time({
	od <- sapply( samplesInCohort , function(u){
		if (!(u%in% samplesInCohort_donors)) return(0)
		i <- which(WW$donor==u)
		sum(WW$infectorProbability[i] ) 
	})
	#})
	
	age <- sapply( sampleDemes[ samplesInCohort ], deme2age )
	stage <- sapply( sampleDemes[ samplesInCohort ], deme2stage )
	
	amat = agemat <- matrix( 0, nrow=4, ncol=4)
	for (k in 1:length(WW$donor)){
		u <- WW$donor[k]
		v <- WW$recip[k]
		au <- deme2age( sampleDemes[ u] )
		av <- deme2age( sampleDemes[v] )
		amat[ au, av ] <- amat[ au, av ] + WW$infectorProbability[k]
	}
	rownames(amat) = colnames(amat) <- 1:4
	
	#load distance mat
	.dfn <-  strsplit( fn, '/' )[[1]]
	dfn <- paste( sep='/', paste(sep='-', .dfn[1], 'distances'), .dfn[2])
	load(dfn) #D
	DD <- D
	D1 <- daytree$tip.label[ D[1,] ]
	D2 <- daytree$tip.label[ D[2,] ]
	DD <- DD[ , (D1 %in% samplesInCohort) & (D2 %in% samplesInCohort) ]
	nbrhoodSize <- setNames( rep(0, length(samplesInCohort)), samplesInCohort )
	agemat2 <- matrix( 0, nrow = 4, ncol = 4 )
	for (k in 1:ncol(DD)){
		if (DD[3,k] < CL_THRESHOLD ) {
			u <- daytree$tip.label[ DD[1,k] ]
			v <- daytree$tip.label[ DD[2,k] ]
			au <- deme2age( sampleDemes[u] )
			av <- deme2age( sampleDemes[v])
			agemat2[au, av] <- agemat2[ au,av] + 1
			nbrhoodSize[u] <- nbrhoodSize[u] + 1
		}
		
	}
	print(date())
	print(fn)
	
	list(agemat = amat
	 , agemat2 = agemat2
	 , tab = data.frame(od=od, age=age, stage=stage , nbrhoodSize = nbrhoodSize)
	 )
}

if (F)
{
fns2agemat_tables <- lapply( ratesBaselineFNS, fn2ageMatrix_table ) 
save(fns2agemat_tables, file = 'age-anal0.0.RData')
} else{
	load('age-anal0.0.RData')
}

mat2trprob <- function(mat) {mat  / rowSums(mat) }
mat2assortmat <- function(mat){
	rs <- rowSums(mat)
	cs <- colSums(mat)
	s <- sum(mat)
	M <- matrix(0, nrow=nrow(mat), ncol = ncol(mat))
	A <- matrix(0, nrow=nrow(mat), ncol = ncol(mat))
	k <- nrow(mat)
	for (i in 1:k) for (j in 1:k){
		M[i,j] <- rs[i] * cs[j] / s
		A[i,j] <- ( mat[i,j] - M[i,j] )/ M[i,j]
	}
	A
}

mat2assortCoef <- function(mat)
{
	mat <- mat / sum(mat )
	rs <- rowSums(mat)
	cs <- colSums(mat)
	sum(diag(mat) - rs*cs) / (1 - sum(rs*cs))
}

if (T)
{
	# check 'true' assort: 
	o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
	tfgy <- .tfgy( o )
	FF <- tfgy[[2]][[length(times_day)]]
	fmat <- matrix(0, nrow = 4, ncol = 4)
	for (k in 1:120) for (l in 1:120){
		fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] <- fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] + FF[k,l]
	}
	BASELINE_ASSRTCOEF <- mat2assortCoef( fmat )

	# make aggregate mixing matrices across all sims
	ag_amat <- fns2agemat_tables[[1]]$agemat
	for (i in 2:length(fns2agemat_tables)){
		ag_amat <- ag_amat + fns2agemat_tables[[i]]$agemat
	}
	# repeat this for clustering
	ag_amat2 <- fns2agemat_tables[[1]]$agemat2
	for (i in 2:length(fns2agemat_tables)){
		ag_amat <- ag_amat + fns2agemat_tables[[i]]$agemat2
	}
	# plot these with assort values also 
	ag_assrt_amat <- mat2assortmat( ag_amat )
	ag_assrt_amat2 <- mat2assortmat( ag_amat2 )
	X11()
	levelplot( ag_amat , main = 'SA mixing matrix')
	X11()
	levelplot( ag_assrt_amat , main = 'SA assortativity')
	X11()
	levelplot( ag_amat2 , main = 'Clustering mixing matrix')
	X11()
	levelplot( ag_assrt_amat2, main = 'Clustering assortativity')

	# compute assort coef for each sim; plot results with 'truth'
	assrt_coefs <- sapply( fns2agemat_tables, function(o){
		mat2assortCoef( o$agemat )
	})
	# repeat this for clustering 
	assrt_coefs2 <- sapply( fns2agemat_tables, function(o){
		mat2assortCoef( o$agemat2 )
	})
	X11()
	boxplot( list( assrt_coefs, assrt_coefs2 ), ylim = c(0.01, 1.25*BASELINE_ASSRTCOEF)
	 , main='Estimated assortativity (1=SA,2=clustering)' 
	 , log = 'y')
	abline( h = BASELINE_ASSRTCOEF, col = 'red')
	#~ p + scale_x_discrete(breaks=c("0.5","1","2"),
	#~         labels=c("Dose 0.5", "Dose 1", "Dose 2"))	
}

if (F)
{
	# for each sims, od ~ age + stage; Q: does age have an indep effect controlling for stage? 
	# A: no, disappears if model includes itxn terms of form stageXage; age does have itxn effect; 
	for (i in 1:length(fns2agemat_tables)){
#~ 	for (i in 1:1){
		tab <- fns2agemat_tables[[i]]$tab
		#lm( od ~ age, data = tab )
		#summary( lm( scale(od) ~ as.factor(age), data = tab ) )
		#summary( lm( scale(log(1+od)) ~ as.factor(age), data = tab ) )
		#require(MASS)
		#boxcox( od ~ as.factor(age) )
		#for ( a in 1:4){
		#	print(summary( tab$od[tab$age==a]))
		#}
		#summary( lm( scale(od) ~ as.factor(age) + as.factor(stage), data = tab ) )
		
		print(
		summary( lm( scale(od) ~ as.factor(age) + as.factor(stage) + as.factor(stage)*as.factor(age), data = tab ) )
		)
	}
	
	# repeat for clustering 
	for (i in 1:length(fns2agemat_tables)){
#~ 	for (i in 1:1){
		tab <- fns2agemat_tables[[i]]$tab
		#lm( od ~ age, data = tab )
		#summary( lm( scale(od) ~ as.factor(age), data = tab ) )
		#summary( lm( scale(log(1+od)) ~ as.factor(age), data = tab ) )
		#require(MASS)
		#boxcox( od ~ as.factor(age) )
		#for ( a in 1:4){
		#	print(summary( tab$od[tab$age==a]))
		#}
		#summary( lm( scale(od) ~ as.factor(age) + as.factor(stage), data = tab ) )
		
		print(
		summary( lm( scale(nbrhoodSize) ~ as.factor(age) + as.factor(stage) + as.factor(stage)*as.factor(age), data = tab ) )
		)
	}
}

# summarise R2 and proportion of sims that show indep sign effect of age on OD or nbrhood size 
if (T)
{
	.summary2age_pvals <- function(s){
		coefs <- s$coefficients
		coefs[ (grepl('(age)', rownames(coefs), fixed=TRUE) & (!grepl('(stage)', rownames(coefs), fixed=TRUE))) , ncol(coefs ) ] 
	}
	nbrhood_R2 <- c()
	nbrhood_age_ps <- c()
	for (i in 1:length(fns2agemat_tables)){
		tab <- fns2agemat_tables[[i]]$tab
		s <- summary( lm( scale(nbrhoodSize) ~ as.factor(age) + as.factor(stage) + as.factor(stage)*as.factor(age), data = tab ) )
		nbrhood_R2 <- c( nbrhood_R2 , s$r.squared )
		nbrhood_age_ps <- c( nbrhood_age_ps, .summary2age_pvals( s ) )
	}
	
	sa_R2 <- c()
	sa_age_ps <- c()
	for (i in 1:length(fns2agemat_tables)){
		tab <- fns2agemat_tables[[i]]$tab
		s <- summary( lm( scale(od) ~ as.factor(age) + as.factor(stage) + as.factor(stage)*as.factor(age), data = tab ) )
		sa_R2 <- c( sa_R2 , s$r.squared )
		sa_age_ps <- c( sa_age_ps, .summary2age_pvals( s ) )
	}
	
	X11();
	boxplot( sa_R2, nbrhood_R2 , main = 'Comparison of RSS for predicting out degree (left) and cluster size(right)')
	print( 'U test comparing SA RSS and clustering neighborhood size RSS:')
	print((wilcox.test( sa_R2, nbrhood_R2 )) )
	print('SA out degree R2:')
	print(summary( sa_R2))
	print('Clustering neighborhood size R2:')
	print(summary( nbrhood_R2 ))
	print ('proportion of p values for SA out degree that reach significance')
	print( mean ( sa_R2  < .05 ))
	print ('proportion of p values for clustering neighborhood size that reach significance')
	print( mean ( nbrhood_R2  < .05 ))
}

if (F)
{
# check sample prop by age group
	table(fns2agemat_tables[[1]]$tab$age )
	YY <- tfgy[[4]][[length(times_day)]]
	ylab2age <- sapply(DEMES[1:120] , deme2age )
	yByAge <- sapply( 1:4, function(a) sum(YY[which(ylab2age==a)] ) )
	cbind( yByAge, table(fns2agemat_tables[[1]]$tab$age ) )
	( yByAge/ table(fns2agemat_tables[[1]]$tab$age ) )
}

if (F)
{
	
	# check 'true' assort: 
	o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
	tfgy <- .tfgy( o )
	FF <- tfgy[[2]][[length(times_day)]]
	fmat <- matrix(0, nrow = 4, ncol = 4)
	for (k in 1:120) for (l in 1:120){
		fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] <- fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] + FF[k,l]
	}
	
	f_trp <- mat2trprob( fmat )
	f_asmat <- mat2assortmat ( fmat )
	
	trp <- mat2trprob( amat )
	asmat <- mat2assortmat ( amat )
	
	trp2 <- mat2trprob( agemat2 )
	asmat2 <- mat2assortmat ( agemat2 )
	
	require(lattice)
	levelplot( asmat )
	X11(); levelplot( f_asmat )
	
	X11(); levelplot( asmat )
	X11(); levelplot( asmat2 )
}
