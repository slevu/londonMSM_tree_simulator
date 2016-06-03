# rm(list=ls())
getwd()
##--- Analyze assortativity matrices computed on HPC 
fn.mat <- list.files('RData', full.names = TRUE, path = "data/simulations2/age/")
load(fn.mat[1])
str(ll)

##- aggregate matrices
##- matrix infector prob
ag_amat_sa <- ll$agemat_sa

for (i in 2:length(fn.mat)){
  load(fn.mat[i])
  ag_amat_sa <- ag_amat_sa + ll$agemat_sa
}

##- aggregate matrix cluster in a list by threshold
str(ll)
##- empty list of 4 matrices
ag_amat_nb <- rep(list(matrix(0,4,4)), length(ll[["agemat_cl"]]))
names(ag_amat_nb) <- names(ll[["agemat_cl"]])
##- loop
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  for (j in 1:length(ll[["agemat_cl"]])){
    ag_amat_nb[[j]] <- ag_amat_nb[[j]] + ll[["agemat_cl"]][[j]][["agemat2"]]
  }
}

##- Newman's assortativity coefficient
mat2assortCoef <- function(mat)
{
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

assor_mat_sa <- mat2assortmat(ag_amat_sa)
assor_mat_nb <- vapply(ag_amat_nb, mat2assortmat, vector(mode = "list", length = 4))

sum(assor_mat_sa)
#######################-------------------------------- laaaaaaaaa
require(lattice)
# lattice.options(default.theme = standard.theme(color=F))
levelplot( assor_mat, col.regions = heat.colors)
levelplot( assor_mat2, col.regions = heat.colors)

