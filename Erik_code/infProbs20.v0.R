# recomputing inf probs using RITA where available
require(phydynR)

MH <- 20 
newinf <- 2500 # c(1660, 4780)
plwhiv <- 43150 / 2 # c(43510 / 2, 43510 / 1.5)


cd4s.df <- read.csv('ukdrd_data/cd4s.csv' )
res.df <- read.csv('ukdrd_data/resistance.csv')
dem.df <- read.csv('ukdrd_data/demographics.csv')

# drop tips that have long cd4 sample time mismatch 
tre <- read.tree( 'LSD/t100.nwk.result.date.newick')
todrop <- readRDS( 'tips2remove_cd4.rds')
tre2 <- drop.tip( tre, colnames(todrop)[which( todrop[1,]==1 )] )

#~ saveRDS( RITA_valid, file = 'RITASampledCloseToSequence.rds' )
RITA_valid <- readRDS(file = 'RITASampledCloseToSequence.rds' )

res.df2 <- res.df[res.df$testindex %in% tre2$tip.label,]
pids <- res.df2$patientindex[match(tre2$tip.label, res.df2$testindex)] 
dem.df2 <- dem.df[ dem.df$patientindex %in% pids , ]
cd4s <- todrop[ 2, tre2$tip.label] 
rita <- dem.df2$ritarecent[ match( pids , dem.df2$patientindex) ]
ehis <- rita
ehis[ehis==99] <- NA
ehis[ehis==1] <- TRUE
ehis[ehis==0] <- FALSE
ehis[ !RITA_valid$RITASampledCloseToSequence[ match( pids, RITA_valid$patientID)] ] <- NA 

sts <- setNames( sapply( as.character( res.df2$dbsample_my) , function( sdt ){
	as.Date( paste(sep='/', '15', sdt ), format='%d/%m/%Y')
}) , res.df2$testindex )

W <- phylo.source.attribution.hiv( tre2,
	sts/365, # years
	cd4s = cd4s,
	ehi = ehis,
	numberPeopleLivingWithHIV = plwhiv,
	numberNewInfectionsPerYear = newinf,
	maxHeight = MH,
	res = 1e3,
	treeErrorTol = 2
)


#####


# out degree 
sids <- tre2$tip.label
od <- sapply( sids, function(u) sum( W$infectorProbability[ W$donor==u ] ))
stage <- sapply( 1:length(sids), function(k){
	cd4 <- cd4s[k]
	ehi <- ehis[k] 
	if (is.na(ehi) | ehi==0){
		if (is.na(cd4)) return(NA)
		if (cd4 > 500 ) return(2)
		if (cd4 > 350 ) return(3)
		if (cd4 > 200 ) return(4)
		return(5)
	} else{
		return(1)
	}
})
od_by_stage <- lapply( 1:5, function(g) {
	x <- od[stage==g] 
	x <- x[!is.na(x)]
	x <- x[x > 0 ]
})
boxplot( lapply( od_by_stage, function(x) log(1e-6+x) ) )
boxplot(  od_by_stage, log = 'y' )

#~ > summary( lm ( od ~ sqrt( cd4s )  ) )
#~ 
#~ Call:
#~ lm(formula = od ~ sqrt(cd4s))
#~ 
#~ Residuals:
#~      Min       1Q   Median       3Q      Max 
#~ -0.05725 -0.02894 -0.02287 -0.01590  2.03858 
#~ 
#~ Coefficients:
#~              Estimate Std. Error t value Pr(>|t|)    
#~ (Intercept) 0.0007908  0.0059181   0.134    0.894    
#~ sqrt(cd4s)  0.0013198  0.0003089   4.273 1.99e-05 ***
#~ ---
#~ Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#~ 
#~ Residual standard error: 0.1176 on 3348 degrees of freedom
#~   (303 observations deleted due to missingness)
#~ Multiple R-squared:  0.005423,	Adjusted R-squared:  0.005126 
#~ F-statistic: 18.25 on 1 and 3348 DF,  p-value: 1.986e-05

age <- 2012 - dem.df2$dob_y[ match( pids , dem.df2$patientindex ) ] 
summary( lm( scale(od) ~ scale(age)*scale(sqrt(cd4s)) ) )
summary( lm( scale(od) ~ scale(age)+scale(sqrt(cd4s)) ) )
#~ > summary( lm( scale(od) ~ scale(age)+scale(sqrt(cd4s)) ) )
#~ 
#~ Call:
#~ lm(formula = scale(od) ~ scale(age) + scale(sqrt(cd4s)))
#~ 
#~ Residuals:
#~     Min      1Q  Median      3Q     Max 
#~ -0.5409 -0.2654 -0.1843 -0.0958 17.1481 
#~ 
#~ Coefficients:
#~                    Estimate Std. Error t value Pr(>|t|)    
#~ (Intercept)       -0.007153   0.016984  -0.421  0.67365    
#~ scale(age)        -0.092702   0.017612  -5.264  1.5e-07 ***
#~ scale(sqrt(cd4s))  0.051390   0.017465   2.942  0.00328 ** 
#~ ---
#~ Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#~ 
#~ i <- which (od > 0 )
#~ summary( lm( log(od[i]) ~ scale(age[i])+scale(sqrt(cd4s[i])) ) )
#~ 

summary( lm( scale(od) ~ as.factor( stage )  ) )
#~ > summary( lm( scale(od) ~ as.factor( stage )  ) )
#~ 
#~ Call:
#~ lm(formula = scale(od) ~ as.factor(stage))
#~ 
#~ Residuals:
#~     Min      1Q  Median      3Q     Max 
#~ -0.3088 -0.2325 -0.1683 -0.1340 17.1004 
#~ 
#~ Coefficients:
#~                   Estimate Std. Error t value Pr(>|t|)
#~ (Intercept)        -0.2147     0.6983  -0.308    0.758
#~ as.factor(stage)2   0.3085     0.6991   0.441    0.659
#~ as.factor(stage)3   0.2323     0.6993   0.332    0.740
#~ as.factor(stage)4   0.1681     0.6990   0.240    0.810
#~ as.factor(stage)5   0.1338     0.6990   0.191    0.848
#~ 
#~ Residual standard error: 0.9875 on 3345 degrees of freedom
#~   (303 observations deleted due to missingness)
#~ Multiple R-squared:  0.004666,	Adjusted R-squared:  0.003476 
#~ F-statistic: 3.921 on 4 and 3345 DF,  p-value: 0.003526
