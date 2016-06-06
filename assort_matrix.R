# rm(list=ls())
getwd()
##--- Analyze assortativity matrices computed on HPC 
fn.mat <- list.files('RData', full.names = TRUE, path = "data/simulations2/age/")
load(fn.mat[1])
str(ll)

##- aggregate matrices
##- matrix of infector probs
ag_amat_sa <- ll$agemat_sa

for (i in 2:length(fn.mat)){
  load(fn.mat[i])
  ag_amat_sa <- ag_amat_sa + ll$agemat_sa
}

##- aggregate matrix cluster in a list by threshold
str(ll)
##- empty list of 4 matrices
ag_amat_nb <- rep(list(matrix(0,4,4)), length(ll[["agemat_cl"]]))
names(ag_amat_nb) <- names(ll[["agemat_cl"]]) # treshold names
##- loop
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  for (j in 1:length(ll[["agemat_cl"]]) ){
    ag_amat_nb[[j]] <- ag_amat_nb[[j]] + ll[["agemat_cl"]][[j]][["agemat2"]]
  }
}

###--- utility functions
##- Newman's assortativity coefficient
mat2assortCoef <- function(mat){
  mat <- mat / sum(mat )
  rs <- rowSums(mat)
  cs <- colSums(mat)
  sum(diag(mat) - rs*cs) / (1 - sum(rs*cs))
}
##- Assortivity matrix (difference from null expectation under random linking)
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

###--- apply assort matrix
assor_mat_sa <- mat2assortmat(ag_amat_sa)
assor_mat_nb <- lapply(ag_amat_nb, mat2assortmat)

###--- plot
require(lattice)
##- SA
levelplot( assor_mat_sa, col.regions = heat.colors)
##- neighborhood
plots_grid <- function(a){
  ## empty list of plots
  p <-  vector("list", length(a)) 
# par(mfrow = c(length(a)/2, 2), bty = 'n') # useless for lattice
for (i in 1:length(a)){
  p[[i]] <- levelplot( a[[i]], col.regions = heat.colors, main = names(a)[i])
}
  return(p)
}

plots_grid(assor_mat_nb)
##- can't find automatic layout, do it manually:
print(p[[1]], split = c(1,1,2,2), more = TRUE)
print(p[[2]], split = c(2,1,2,2), more = TRUE)
print(p[[3]], split = c(1,2,2,2), more = TRUE)
print(p[[4]], split = c(2,2,2,2), more = FALSE)

###--- apply assort coef
coef_mat_sa <- mat2assortCoef(ag_amat_sa)
coef_mat_nb <- lapply(ag_amat_nb, mat2assortCoef)

##-- true assortativity from model
deme2age <- function(deme){as.numeric( substr(regmatches( deme , regexpr( '\\.age[0-9]', deme )), 5,5) ) }
deme2stage <- function(deme){as.numeric( substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6,6) ) }

require(phydynR)
source('model0.R')
o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
tfgy <- .tfgy( o )
FF <- tfgy[[2]][[length(times_day)]]
fmat <- matrix(0, nrow = 4, ncol = 4)
for (k in 1:120) for (l in 1:120){
  fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] <- fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] + FF[k,l]
}
BASELINE_ASSRTCOEF <- mat2assortCoef( fmat )
print( 'True assortativity in baseline simulation')
print(BASELINE_ASSRTCOEF)
##- "true" heatmap
levelplot( mat2assortmat( fmat ), col.regions = heat.colors)

##- compare true with estimated coefs
# compute assort coef for each sim; plot results with 'truth'


assrt_coefs_sa <- vector(mode="numeric", length = length(fn.mat))

for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  assrt_coefs_sa[i] <- mat2assortCoef(ll$agemat_sa)
}

# repeat this for clustering 
# ##- empty list of 4 matrices
assrt_coefs_nb <- rep(list(vector(mode="numeric", length = length(fn.mat))), length(ll[["agemat_cl"]]))
names(assrt_coefs_nb) <- names(ll[["agemat_cl"]]) # treshold names
##- loop
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  for (j in 1:length(ll[["agemat_cl"]]) ){
    assrt_coefs_nb[[j]][i] <-   mat2assortCoef(ll[["agemat_cl"]][[j]][["agemat2"]])
  }
}

##- transfrom in df
df <- cbind("SA" = assrt_coefs_sa, 
            do.call(cbind.data.frame, assrt_coefs_nb) )
str(df)
##################------------------- llllaaa
temp = reshape(df, direction="long", varying=1:5, sep="")
##- boxplot
boxplot( list( assrt_coefs, assrt_coefs2 ), ylim = c(0.01, 1.25*BASELINE_ASSRTCOEF)
         , main='Estimated assortativity (1=SA,2=clustering)' 
         , log = 'y')
abline( h = BASELINE_ASSRTCOEF, col = 'red')