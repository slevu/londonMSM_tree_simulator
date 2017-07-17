# rm(list=ls())
####---- include ----
detail_knitr <- TRUE

##---- libs ----
library(reshape2)
library("ggplot2")
# theme_set(theme_bw())
library(scales)
library(gridExtra)
require(lattice)

##---- format exponential ----
# options(scipen = -100)
# options(scipen = 0)

##--- compute age matrix for UCSD cluster size  ---
##---- load data ----
source('load_sim_results.R')

# names(cw_Baseline0) ; head(cw_Baseline0[[1]][[1]])
  ##- cluster and nbhsize lists only
  c_Base <- cw_Baseline0[-1]
  ##- load EdgeList
  source("functions.R")

if(F){
##- For ucsd clustering
  ##- empty list of 4 matrices
  ##list(matrix(0,4,4))
  list_agmat_cl <- rep(list(rep(list(matrix(0,4,4)), length = length(c_Base[[1]]) )), length(c_Base))
  names(list_agmat_cl) <- names(c_Base)
  # treshold names
  ##- loop
system.time(
  for (i in 1:length(c_Base)){ # index threshold
    for (j in 1:length(c_Base[[i]]) ){ # index sim
      l <- c_Base[[i]][[j]]
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
      names(list_agmat_cl[[i]])[j] <- names(c_Base[[i]])[j]
    }
  }
) # 720s !
  
  saveRDS(list_agmat_cl, file =  paste(path.results, 'model1-sim_ucsd/list.agmat.clus.rds', sep = '/'))
  
  ##- rm functions and large file
  rm(list = lsf.str(), cw_Baseline0)
  
} else {
  list_agmat_cl <- readRDS(paste(path.results, 'model1-sim_ucsd/list.agmat.clus.rds', sep = '/'))
}
# options(scipen = -100)
names(list_agmat_cl) <- as.character( as.numeric(names(list_agmat_cl) ) )
# options(scipen = 0)
# names(list_agmat_cl); head(names(list_agmat_cl[[1]])); lapply(list_agmat_cl, function(x) x[1])

##--- Analyze age matrices computed on HPC 
fn.mat <- list.files('RData', full.names = TRUE, path = paste(path.results, 'model1_age', sep = '/'))

##---- aggregate ----
##- aggregate matrix cluster in a list by threshold
##- ucsd cluster
ag_amat_cl <- lapply(list_agmat_cl, function(x) Reduce('+', x))

##-  age matrix of infector probs
load(fn.mat[1]) 
# str(ll)
ag_amat_sa <- ll$agemat_sa
for (i in 2:length(fn.mat)){
  load(fn.mat[i])
  ag_amat_sa <- ag_amat_sa + ll$agemat_sa
}

##- nbrhood size
ag_amat_nb <- rep(list(matrix(0,4,4)), length(ll[["agemat_cl"]]))
# options(scipen = -100)
names(ag_amat_nb) <- as.character( as.numeric(names(ll[["agemat_cl"]])) ) # treshold names
# options(scipen = 0)
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  for (j in 1:length(ll[["agemat_cl"]]) ){
    ag_amat_nb[[j]] <- ag_amat_nb[[j]] + ll[["agemat_cl"]][[j]][["agemat2"]]
  }
}

##---- utility functions ----
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

##---- apply assort matrix ----
assor_mat_sa <- mat2assortmat(ag_amat_sa)
assor_mat_nb <- lapply(ag_amat_nb, mat2assortmat)
assor_mat_cl <- lapply(ag_amat_cl, mat2assortmat)

##---- apply assort coef globally ----
coef_mat_sa <- mat2assortCoef(ag_amat_sa)
coef_mat_nb <- lapply(ag_amat_nb, mat2assortCoef)
coef_mat_cl <- lapply(ag_amat_cl, mat2assortCoef)

##---- normalize ----
normalit <- function(m){
  M <- (m - min(m))#/(max(m)-min(m))
  dimnames(M) <- list(1:nrow(m), 1:ncol(m))
  return(M)
}
norm_assor_mat <- lapply(c(sa = list(assor_mat_sa), cl = assor_mat_cl, nb = assor_mat_nb), normalit)
coefs <- c(sa = coef_mat_sa, cl = unlist(coef_mat_cl), nb = unlist(coef_mat_nb))
mx <-  max(sapply(norm_assor_mat, function(m) max(m))) # max of all matrices
bkpoints <- seq(0, mx, length.out=20) # same breakpoints for all levelplots
TOTO <- do.call(rbind, lapply(1:length(norm_assor_mat), function(x){
  df <- as.data.frame(as.table(norm_assor_mat[[x]]))
  df$group <- paste0(letters[x], ' (r =', round(coefs[x],2), ')')
  return(df)
} ))
levelplot(Freq ~ Var1 * Var2 | group, data = TOTO, col.regions = heat.colors,
          par.settings = list(strip.background=list(col="white")),
          xlab = '', ylab = '', as.table = TRUE ) # order with as.table or index.cond=list()

##---- heatmap 1 ----
##- SA
# print("Assortativity of infector probs by age")
p_sa <- levelplot( norm_assor_mat[['sa']], col.regions = heat.colors,
           main = 'SA', 
           sub = paste('r = ', round(coef_mat_sa,2)), 
           xlab = '', ylab = '', at = bkpoints )

# p_sa

##---- heatmap 2 ----
##- neighborhood
plots_grid <- function(mat, coef){
  ## empty list of plots
  p <-  vector("list", length(mat)) 
  # par(mfrow = c(length(a)/2, 2), bty = 'n') # useless for lattice
  for (i in 1:length(mat)){
    p[[i]] <- levelplot( mat[[i]], col.regions = heat.colors, main = names(mat)[i], sub = paste('r = ', round(coef[[i]],2)), xlab = '', ylab = '', at = bkpoints )
  }
  return(p)
}

p_nb <- plots_grid(norm_assor_mat[grep('nb', names(norm_assor_mat))], coef_mat_nb[1:4])
# print("Assortativity of neighborhood size by age")
##- plot
# do.call(grid.arrange, c(p_nb, ncol = ceiling(length(p_nb)/2)) )

##---- heatmap 3 ----
##- ucsd
p_cl <- plots_grid(norm_assor_mat[grep('cl', names(norm_assor_mat))], coef_mat_cl[1:4])
# print("Assortativity of cluster size by age")
# do.call(grid.arrange, c(p_cl, ncol = ceiling(length(p_cl)/2)) )

##- arrange all plots
all_p <- c(list(p_sa), p_cl, p_nb)
##- rename main title with a,b,c and r newman, small text plain, par.settings
theme1 <-
list(layout.heights =
       list(top.padding = 1,
            main.key.padding = .5,
            key.axis.padding = 0,
            axis.xlab.padding = 0,
            xlab.key.padding = 0,
            key.sub.padding = .5,
            bottom.padding = 1),
     layout.widths =
       list(left.padding = .5,
            key.ylab.padding = 0,
            ylab.axis.padding = 0,
            axis.key.padding = 0,
            right.padding = .5))

all_p <- lapply(1:length(all_p), function(x){
  all_p[[x]]$main <- list(label = paste0('(', letters[x], ')'), cex = 1)
  all_p[[x]]$sub <- list(label = all_p[[x]]$sub, cex = 1, font = 1)
  all_p[[x]]$par.settings <- theme1
  return( all_p[[x]])
}   )

do.call(grid.arrange, c( all_p, ncol = 3 ))


##---- true assortativity from model0 ----
deme2age <- function(deme){as.numeric( substr(regmatches( deme , regexpr( '\\.age[0-9]', deme )), 5,5) ) }
deme2stage <- function(deme){as.numeric( substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6,6) ) }

##- run model
if(FALSE){
  require(phydynR)
source('HPC_code/model1.R')
o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
tfgy <- .tfgy( o )
FF <- tfgy[[2]][[length(times_day)]]
fmat <- matrix(0, nrow = 4, ncol = 4)
for (k in 1:120) for (l in 1:120){
  fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] <- fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] + FF[k,l]
}
saveRDS(fmat, file = paste(path.results, 'model1-sim_ucsd/model1.true_agemat.rds', sep = '/'))
} else {
  fmat <- readRDS(file = paste(path.results, 'model1-sim_ucsd/model1.true_agemat.rds', sep = '/'))
}
##- True assortativity in baseline simulation
BASELINE_ASSRTCOEF <- mat2assortCoef( fmat )
print(BASELINE_ASSRTCOEF)
##- "true" heatmap
# levelplot( mat2assortmat( fmat ), col.regions = heat.colors)


##---- estimated coefs ----
# compute assort coef for each sim; plot results with 'truth'
assrt_coefs_sa <- vector(mode="numeric", length = length(fn.mat))
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  assrt_coefs_sa[i] <- mat2assortCoef(ll$agemat_sa)
}

##- For neighborhood
##- empty list of 4 matrices
assrt_coefs_nb <- rep(list(vector(mode="numeric", length = length(fn.mat))), length(ll[["agemat_cl"]]))
# options(scipen = -100)
names(assrt_coefs_nb) <- as.character(as.numeric ( names(ll[["agemat_cl"]]) ) ) # treshold names
# options(scipen = 0)

##- loop
for (i in 1:length(fn.mat)){
  load(fn.mat[i])
  for (j in 1:length(ll[["agemat_cl"]]) ){
    assrt_coefs_nb[[j]][i] <-   mat2assortCoef(ll[["agemat_cl"]][[j]][["agemat2"]])
  }
}

##- For ucsd
##- empty list of 4 matrices
assrt_coefs_cl <-  lapply(list_agmat_cl, function(x) sapply(x, mat2assortCoef))
##---- stop ----

##---- boxplot of assort coef ----


##- transform in df
df <- cbind("SA" = assrt_coefs_sa, 
            "NB" = do.call(cbind.data.frame, assrt_coefs_nb),
            "CL" = do.call(cbind.data.frame, assrt_coefs_cl) )
# str(df)
## with reshape2
a <- melt(df)
## extra columns for facets
a$method <- substr(a$variable, 1, 2)
##- add '' for threshold of SA method
a$thr <- c(rep('NA', length(a[a$method == 'SA', 'method'])), regmatches(a$variable, regexpr("\\d\\..*|\\d*e[-+]?.*",  a$variable)) )
# table(a$thr); table(a$method)
# str(a)
##- keep only 4 thersholds
a <-a[!(a$thr %in% c('1e-05','1e-04')),]
##- force scientific notation
fancy_scientific <- function(l) {
  if( is.na(as.numeric(l)) ) parse(text = l) else {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
  }
}

##---- boxplot 1 ----
library(cowplot)
bp1 <- ggplot(a, aes(reorder(thr, as.numeric(thr)), value, color = method)) + geom_boxplot()
bp2 <- bp1  + 
  facet_grid(~ method, scales = "free", space = "free", labeller=labeller(method = c(SA = "SA", CL = "Cluster", NB = "Neighborhood")))  + 
  geom_hline(aes(yintercept = BASELINE_ASSRTCOEF, colour = "true coefficient"), linetype = "dashed") + 
  xlab("Distance threshold") + ylab("Assortativity coefficient")

  # theme(legend.position="top", legend.title=element_blank(), strip.text.x = element_text(size = 12)) + guides(fill=FALSE) +
bp3 <- bp2 + theme(legend.position="none") + 
  background_grid() + # + theme_bw()
theme(strip.background = element_blank())
bp3
# bp2 %+% a[a$thr != '1e-05',] + scale_x_discrete(labels = fancy_scientific) + theme(legend.position="none") + 
# background_grid() + theme(strip.background = element_blank())

##---- boxplot 2 ----
## With scale that adapts to transformation
bp3 %+% a[a$value >= 0,] + coord_trans(y = "log")

##---- boxplot age base ----
# str(a)
# str(df)
##- to keep same range
lim <- range(a$value, 1.1*BASELINE_ASSRTCOEF)

# par(mfrow=c(1,3) ) #, bty = 'n')   
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), 
       widths=c(1,1.5,1.5), heights=c(1,1,1))

boxplot(a[a$method == 'SA',]$value ~ a[a$method == 'SA',]$thr, ylim = lim)
title(main = 'SA', font.main = 1, ylab = 'Assortativity coefficient', cex.lab = 1.2)
abline(h = BASELINE_ASSRTCOEF, lty = 3)

boxplot(a[a$method == 'CL',]$value ~ a[a$method == 'CL',]$thr, ylim = lim, yaxt="n")
title(main = 'Cluster', xlab = 'Distance threshold', cex.lab = 1.2,  font.main = 1)
abline(h = BASELINE_ASSRTCOEF, lty = 3)

boxplot(a[a$method == 'NB',]$value ~ a[a$method == 'NB',]$thr, ylim = lim, yaxt="n")
title(main = 'Neighborhood', xlab = 'Distance threshold', cex.lab = 1.2,  font.main = 1)
abline(h = BASELINE_ASSRTCOEF, lty = 3)

# dev.off()
# ?abline; ?layout; ?boxplot

