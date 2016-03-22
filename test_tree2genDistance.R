### test tree2genDistance.R

## simulate a dated tree as in PhyDynR example
require(phydynR)
dm.det <- build.demographic.process(births = c(I = 'parms$beta * I' )
                                    , deaths = c(I = 'parms$gamma * I' )
                                    , parameterNames=c('beta', 'gamma') 
                                    , rcpp= FALSE
                                    , sde = FALSE)

show.demographic.process( dm.det
                          , theta = list( beta = 1.5, gamma = 1 )
                          , x0  = c( I = 1 )
                          , t0 = 0
                          , t1 = 10 
) 

## simulate tree
tree <- sim.co.tree(   list( beta = 1.5, gamma = 1 )
, dm.det
, x0  = c(I = 1 )
, t0 = 0
, sampleTimes = seq(10, 15, length.out = 10)
, res = 1000
) 

names(tree)
plot(tree)

source("tree2genDistance.R")
# (year * susbt/site/year * n site) = number of substitutions 
D <-  .tree2genDistance(tree, mu = .0015, seqlength = 1e3, MH=20 )

