# rm(list=ls())

##---- compute age matrix for UCSD cluster size  ----
##---- load data ----
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
##- load EdgeList
source("functions.R")

##- For ucsd clustering
if(FALSE){
  ##- empty list of 4 matrices
  ##list(matrix(0,4,4))
  list_agmat_cl <- rep(list(rep(list(matrix(0,4,4)), length = length(cw_Baseline0[[1]]))), length(cw_Baseline0))
  names(list_agmat_cl) <- names(cw_Baseline0)
  # treshold names
  ##- loop
  for (i in 1:length(cw_Baseline0)){ # index threshold
    for (j in 1:length(cw_Baseline0[[i]]) ){ # index sim
      l <- cw_Baseline0[[i]][[j]]
      ##- create vector of age (name) assigned to cluster number
      x <- l$ClusterID
      names(x) <- l$age
      ##- compute edge list
      el <- EdgeList(x)
      ##- empty matrix
      agemat <- matrix( 0, nrow = length(unique(el[,'from'])), ncol = length(unique(el[,'to'])) )
      ##- populate matrix
      for (k in 1:dim(el)[1]){
        from <- el[k, 1]
        to <- el[k, 2]
        agemat[from, to] <- agemat[from, to] + 1
      }
      list_agmat_cl[[i]][[j]] <-  agemat
      names(list_agmat_cl[[i]])[j] <- names(cw_Baseline0[[i]])[j]
    }
  }
  # saveRDS(list_agmat_cl, file = "data/sim_ucsd_results2/list.agmat.clus.rds" )
} else {
  list_agmat_cl <- readRDS("data/sim_ucsd_results2/list.agmat.clus.rds")
}
# str(list_agmat_cl[[1]])

##- rm fns and large file
rm(list = lsf.str(), cw_Baseline0)

##- aggregate matrix cluster in a list by threshold
##- empty list of 4 matrices
ag_amat_cl <- rep(list(matrix(0,4,4)), length(list_agmat_cl))
names(ag_amat_cl) <- names(list_agmat_cl) # treshold names
##- loop
for (i in 1:length(list_agmat_cl) ){
  for (j in 1:length(list_agmat_cl[[i]]) ){
    ag_amat_cl[[i]] <- ag_amat_cl[[i]] + list_agmat_cl[[i]][[j]]
  }
}


##--- Analyze assortativity matrices computed on HPC 
fn.mat <- list.files('RData', full.names = TRUE, path = "data/simulations2/age/")
##- load first list
load(fn.mat[1]) 
# str(ll)

##- aggregate matrices
##- matrix of infector probs
ag_amat_sa <- ll$agemat_sa

for (i in 2:length(fn.mat)){
  load(fn.mat[i])
  ag_amat_sa <- ag_amat_sa + ll$agemat_sa
}

##- aggregate matrix cluster in a list by threshold
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
assor_mat_cl <- lapply(ag_amat_cl, mat2assortmat)

###--- plots ---
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

p <- plots_grid(assor_mat_nb)
##- can't find automatic layout, do it manually:
print(p[[1]], split = c(1,1,2,2), more = TRUE)
print(p[[2]], split = c(2,1,2,2), more = TRUE)
print(p[[3]], split = c(1,2,2,2), more = TRUE)
print(p[[4]], split = c(2,2,2,2), more = FALSE)

##- ucsd
p <- plots_grid(assor_mat_cl)
##- can't find automatic layout, do it manually:
print(p[[1]], split = c(1,1,2,2), more = TRUE)
print(p[[2]], split = c(2,1,2,2), more = TRUE)
print(p[[3]], split = c(1,2,2,2), more = TRUE)
print(p[[4]], split = c(2,2,2,2), more = FALSE)

###--- apply assort coef
coef_mat_sa <- mat2assortCoef(ag_amat_sa)
coef_mat_nb <- lapply(ag_amat_nb, mat2assortCoef)
coef_mat_cl <- lapply(ag_amat_cl, mat2assortCoef)

##-- true assortativity from model
deme2age <- function(deme){as.numeric( substr(regmatches( deme , regexpr( '\\.age[0-9]', deme )), 5,5) ) }
deme2stage <- function(deme){as.numeric( substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6,6) ) }

##- run model
if(FALSE){
  require(phydynR)
source('model0.R')
o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
tfgy <- .tfgy( o )
FF <- tfgy[[2]][[length(times_day)]]
fmat <- matrix(0, nrow = 4, ncol = 4)
for (k in 1:120) for (l in 1:120){
  fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] <- fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] + FF[k,l]
}
} else {
  # saveRDS(fmat, file = "data/sim_ucsd_results2/model0.agemat.rds")
  fmat <- readRDS(file = "data/sim_ucsd_results2/model0.agemat.rds")
}
##- True assortativity in baseline simulation
BASELINE_ASSRTCOEF <- mat2assortCoef( fmat )
print(BASELINE_ASSRTCOEF)
##- "true" heatmap
levelplot( mat2assortmat( fmat ), col.regions = heat.colors)

##--- compare true with estimated coefs ---
# compute assort coef for each sim; plot results with 'truth'
assrt_coefs_sa <- vector(mode="numeric", length = length(fn.mat))
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  assrt_coefs_sa[i] <- mat2assortCoef(ll$agemat_sa)
}

##- For neighborhood
##- empty list of 4 matrices
assrt_coefs_nb <- rep(list(vector(mode="numeric", length = length(fn.mat))), length(ll[["agemat_cl"]]))
names(assrt_coefs_nb) <- names(ll[["agemat_cl"]]) # treshold names
##- loop
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  for (j in 1:length(ll[["agemat_cl"]]) ){
    assrt_coefs_nb[[j]][i] <-   mat2assortCoef(ll[["agemat_cl"]][[j]][["agemat2"]])
  }
}

##- For ucsd
##- empty list of 4 matrices
assrt_coefs_cl <- rep(list(vector(mode="numeric", length = length(list_agmat_cl[[1]]))), length(list_agmat_cl))
names(assrt_coefs_cl) <- names(list_agmat_cl) # treshold names
##- loop
for (i in 1:length(list_agmat_cl)){
  for (j in 1:length(list_agmat_cl[[1]]) ){
    assrt_coefs_cl[[i]][[j]] <-   mat2assortCoef(list_agmat_cl[[i]][[j]])
  }
}
## or
# sapply(list_agmat_cl, function(x){
#   lapply(x, mat2assortCoef)
#   })

##-- boxplot of assort coef
##- transform in df
names(assrt_coefs_nb) <- paste0("nb", names(assrt_coefs_nb))
names(assrt_coefs_cl) <- paste0("cl", names(assrt_coefs_cl))


df <- cbind("SA" = assrt_coefs_sa, 
            do.call(cbind.data.frame, assrt_coefs_nb),
            do.call(cbind.data.frame, assrt_coefs_cl))
# str(df)
library(reshape2)
a <- melt(df)

##- boxplot
# boxplot( a$value ~ a$variable, ylim = c(min(a$value), 1.25*BASELINE_ASSRTCOEF)
#          , main='Estimated assortativity (SA,nb=neighborhood, cl=clustering)' 
#          , log = 'y')
# abline( h = BASELINE_ASSRTCOEF, col = 'red')

library(ggplot2)
bp <- ggplot(a, aes(variable, value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = BASELINE_ASSRTCOEF, colour = "true coef"))
bp + theme_bw()


require(lattice)
# lattice.options(default.theme = standard.theme(color=F))
levelplot( mat2assortmat( agemat ), col.regions = heat.colors)

# ##- essai neighborhood size
# cor(l$size, l$nbhsize)
# cor(l$outdegree, l$nbhsize, use = "complete" )
# plot(l$outdegree, l$nbhsize)
# boxplot(l$outdegree ~ l$risk) 
# boxplot(l$nbhsize ~ l$risk) 
# boxplot(l$size ~ l$risk) 