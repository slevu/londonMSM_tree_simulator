# simulate a tree: 
source('model0.R') 

# find the stage for each patient id (type=character); NOTE finding the corresponding element should reference tree$sampleStates, not the input sample states
pid2stage <- function(pid)
{
	deme <- colnames(tree$sampleStates)[which.max(tree$sampleStates[pid,]) ]
	as.numeric( regmatches( deme, regexec( "stage([0-9])", deme) )[[1]][2] )
}

# test 1, time distance, dep variable = number samples with time dist < threshold
D <- cophenetic.phylo( tree )
diag(D) <- Inf
time_threshold <- 10*365 
if (T)
{
	
	odegree <- sapply( rownames(D), function(tl) sum( D[tl,] < time_threshold ))
	odegree <- sapply( 1:nrow(D), function(k) sum( D[k,] < time_threshold ))
	hist( odegree )
	#save.image( file = 'm0.testclusteringAssoc.RData' ) # Too slow
	save( tree, odegree, time_threshold, file = 'm0.testclusteringAssoc.RData'  )
	
	pids <- rownames(D)
	stages <- sapply( pids, pid2stage )
	#lm ( stages ~ odegree  )
#~ 	> summary( lm ( stages ~ odegree  ) )
#~ 
#~ Call:
#~ lm(formula = stages ~ odegree)
#~ 
#~ Residuals:
#~      Min       1Q   Median       3Q      Max 
#~ -2.36310 -1.07332 -0.07332  0.92668  3.08578 
#~ 
#~ Coefficients:
#~              Estimate Std. Error t value Pr(>|t|)    
#~ (Intercept)  3.363099   0.018504  181.75   <2e-16 ***
#~ odegree     -0.048296   0.002369  -20.39   <2e-16 ***
#~ ---
#~ Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#~ 
#~ Residual standard error: 1.256 on 12162 degrees of freedom
#~ Multiple R-squared:  0.03304,	Adjusted R-squared:  0.03296 
#~ F-statistic: 415.6 on 1 and 12162 DF,  p-value: < 2.2e-16
#~ 

#~ > summary( glm ( (odegree>0) ~ stages, family=binomial(link='logit') ) )
#~ 
#~ Call:
#~ glm(formula = (odegree > 0) ~ stages, family = binomial(link = "logit"))
#~ 
#~ Deviance Residuals: 
#~     Min       1Q   Median       3Q      Max  
#~ -2.7426   0.2675   0.3292   0.4043   0.4949  
#~ 
#~ Coefficients:
#~             Estimate Std. Error z value Pr(>|z|)    
#~ (Intercept)  4.16199    0.12280   33.89   <2e-16 ***
#~ stages      -0.42473    0.03203  -13.26   <2e-16 ***
#~ ---
#~ Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#~ 
#~ (Dispersion parameter for binomial family taken to be 1)
#~ 
#~     Null deviance: 5588.4  on 12163  degrees of freedom
#~ Residual deviance: 5398.8  on 12162  degrees of freedom
#~ AIC: 5402.8
#~ 
#~ Number of Fisher Scoring iterations: 6
#~ 


}
