# rm(list=ls())

##---- libs ----
library(reshape2)
library("ggplot2")
theme_set(theme_bw())
library(scales)

##---- compute age matrix for UCSD cluster size  ----
##---- load data ----
if(FALSE){
  cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
# names(cw_Baseline0) ; head(cw_Baseline0[[1]][[1]])
  ##- load EdgeList
  source("functions.R")

##- For ucsd clustering
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
  
  ##- rm functions and large file
  rm(list = lsf.str(), cw_Baseline0)
  
} else {
  list_agmat_cl <- readRDS("data/sim_ucsd_results2/list.agmat.clus.rds")
}
# names(list_agmat_cl); head(names(list_agmat_cl[[1]]))

##--- Analyze age matrices computed on HPC 
fn.mat <- list.files('RData', full.names = TRUE, path = "data/simulations2/age/")

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
names(ag_amat_nb) <- names(ll[["agemat_cl"]]) # treshold names
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

##---- heatmaps ----
require(lattice)
##- SA
print("Assortativity of infector probs by age")
levelplot( assor_mat_sa,
           col.regions = heat.colors)

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
print("Assortativity of neighborhood size by age")
print(p[[1]], split = c(1,1,2,2), more = TRUE)
print(p[[2]], split = c(2,1,2,2), more = TRUE)
print(p[[3]], split = c(1,2,2,2), more = TRUE)
print(p[[4]], split = c(2,2,2,2), more = FALSE)

##- ucsd
p <- plots_grid(assor_mat_cl)
##- can't find automatic layout, do it manually:
print("Assortativity of cluster size by age")
print(p[[1]], split = c(1,1,2,2), more = TRUE)
print(p[[2]], split = c(2,1,2,2), more = TRUE)
print(p[[3]], split = c(1,2,2,2), more = TRUE)
print(p[[4]], split = c(2,2,2,2), more = FALSE)


##---- apply assort coef ----
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
assrt_coefs_cl <-  lapply(list_agmat_cl, function(x) sapply(x, mat2assortCoef))
##---- stop ----

##---- boxplot of assort coef ----
##- transform in df
names(assrt_coefs_nb) <- paste0("nb", names(assrt_coefs_nb))
names(assrt_coefs_cl) <- paste0("cl", names(assrt_coefs_cl))

df <- cbind("SA" = assrt_coefs_sa, 
            "NB" = do.call(cbind.data.frame, assrt_coefs_nb),
            "CL" = do.call(cbind.data.frame, assrt_coefs_cl) )
# str(df)
## with reshape2
a <- melt(df)
## extra columns for facets
a$method <- substr(a$variable, 1, 2)
##- add '' for threshold of SA method
a$thr <- c(rep('NA', length(a[a$method == 'SA', 'method'])), regmatches(a$variable, regexpr("\\d\\.\\d*",  a$variable)) )

##- boxplot
# boxplot( a$value ~ a$variable, ylim = c(min(a$value), 1.25*BASELINE_ASSRTCOEF), main='Estimated assortativity (SA,nb=neighborhood, cl=clustering)', log = 'y')
# abline( h = BASELINE_ASSRTCOEF, col = 'red')

library(ggplot2)
bp1 <- ggplot(a, aes(thr, value, fill = method)) + geom_boxplot()
bp2 <- bp1  + 
  facet_grid(~ method, scales = "free", space = "free", labeller=labeller(method = c(SA = "SA", CL = "UCSD cluster", NB = "Neighborhood")))  + 
  geom_hline(aes(yintercept = BASELINE_ASSRTCOEF, colour = "true coefficient")) + 
  theme(legend.position="top", legend.title=element_blank(), strip.text.x = element_text(size = 12)) + guides(fill=FALSE) +
  xlab("Distance threshold") + ylab("Assortativity coefficient")
bp2

## With scale that adapt to transformation
bp2 %+% a[a$value >= 0,] + coord_trans(y = "log")

 
##---- associations ----
## age vs sizes, degrees at different thr
rm(list=ls())
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )

## change structure in long table (first few for speed)
.n <- 3 
.m <- sample(1:100, .n)
c <- lapply(cw_Baseline0, 
            function(x) do.call(rbind, x[.m]))
rm(cw_Baseline0)

# head(c[[2]])
# names(c)
# str(c)
# lapply(c, function(x) mean(x$nbhsize))

## winsorize outliers
winsorize <- function (x, fraction=.05)
{
  if(length(fraction) != 1 || fraction < 0 ||
     fraction > 0.5) {
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction, 1-fraction), na.rm = TRUE)
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}

cc <- lapply(c, function(x) {
  x[,c(6, 8, 10)] <- lapply(x[,c(6, 8, 10)], winsorize)
  return(x)
  })
str(cc)
library(reshape2)
df <- melt(c,  measure.vars = c(6, 8, 10))

sa <- melt( cc[[1]][, -c(6, 9:10)],  measure.vars = 7)
# head(sa)
rest <- melt(cc,  measure.vars = c(6, 10))
# head(rest)
tot <- rbind( cbind(sa, "L1" = 'NA'), rest[,-c(7,8)])
# head(tot)

df_win <- melt(cc,  measure.vars = c(6, 8, 10))
# head(df)
# str(df)
# sum(is.na(df$value))
table(df$L1)
identical(c[[1]][, "outdegree"], c[[3]][, "outdegree"]) # TRUE
identical(df[df$L1 == '0.001' & df$variable == 'outdegree', 'value'], 
          df[df$L1 == '0.005' & df$variable == 'outdegree', 'value'])

# df <- df[df$L1 == "0.015",]
tapply(df$value, list(df$L1, df$variable), function(x) mean(x, na.rm = T))

p <- ggplot(df, aes(x = factor(age), y = value, fill = variable)) + geom_boxplot() + theme_bw() + theme(legend.position="none", legend.title = element_blank(), strip.text.x = element_text(size = 12)) 
p2 <- p + facet_wrap(~ L1 + variable, nrow = 4, scales = "free_y")
p2
p2 %+% df_win

p3 <- ggplot(tot, aes(factor(age), value, fill = variable)) + geom_boxplot()
p4 <- p3 + facet_wrap(~ variable + L1, nrow = 1, scales = "free_y")
p4

##- at 0.015 only, winsorized
p %+% df[df$L1 == '0.015',] + facet_wrap(~ variable, nrow = 1, scales = "free_y")
p %+% df_win[df_win$L1 == '0.015',] + facet_wrap(~ variable, nrow = 1, scales = "free_y", labeller=labeller(variable = c(outdegree = "Outdegree", size = "Cluster size", nbhsize = "Neighborhood \n size")))


##---- regressions ----
###--- individual bootstrap regressions ---
###--- start function
###- summarize regression on bootstrap
# x = cw_Baseline0[[1]][[1]]; ls = cw_Baseline0, reg =lm; model =  "scale(size) ~ factor(age)"
# y = 'scale(size)';  xs = c('factor(age)')
reg.sum.bs <- function(ls, reg, y, xs, alpha = 0.05, ...){
  
  model <- paste(y, '~', paste(xs, collapse = ' + '))
  
  ## coef by threshold and by tree
  coef <- lapply(ls, function(x){
    lapply(x , function(x){
      coef(summary(reg(formula = model, data = x, ...)))
    })
  })
  # str(coef[[1]][[1]])
  
  ## pvalue by threshold and by tree
  pvalue <- lapply(coef, function(x){
    sapply(x , function(x){
      identity(x[,4])
    })
  })
  
  ##- number of p-value < 0.05
  sum.signif <- sapply(pvalue, function(x){
    apply(x, 1, function(x) sum(x < alpha) / length(x))
  })
  
  ## number of p-value < 0.05 for first x of xs only and without interaction
  sum.signif_x <- sapply(pvalue, function(x){ 
    sum(x[grepl(xs[1], rownames(x), fixed=TRUE) & !grepl(':', rownames(x), fixed=TRUE),] < alpha) / length(x[grepl(xs[1], rownames(x), fixed=TRUE) & !grepl(':', rownames(x), fixed=TRUE), ])
  })
  
  ## parameter by threshold
  param <-  lapply(coef, function(x){
    sapply(x , function(x){
      identity(x[,1])
    })
  })
  
  ## mean of parameter
  mean.parms <- signif(sapply(param, function(x){
    apply(x, 1, mean)
  }), 2)
  
  ## R square, only for lm()
  if(identical(reg, lm)){
    r2 <- lapply(ls, function(x){
      sapply(x , function(x){
        summary(reg(model, data = x))$r.squared
      })
    })
    ## mean R2
    mean.r2 <- signif(sapply(r2, function(x){
      mean(x)
    }), 3)
    
    return(list("model" = model, "mean parameter" = mean.parms, "signif pvalue" = sum.signif,
                "signif pvalue xs[1]" = sum.signif_x, "mean r.squared" = mean.r2)) 
  } else {
    
    return(list("model" = model, "mean parameter" = mean.parms, "signif pvalue" = sum.signif, "signif pvalue xs[1]" = sum.signif_x))
  }
}
###--- end function 


model3_cl <- "scale(size) ~ factor(age)"
model3_sa <- "scale(outdegree) ~ factor(age)"
model3_nb <- "scale(nbhsize) ~ factor(age)"
model4 <- "scale(size) ~ factor(age) + factor(stage) + factor(age)*factor(stage)"

xs1 = c('factor(age)', 'factor(stage)', 'factor(age)*factor(stage)')
## example
## age
a <- reg.sum.bs(ls = cw_Baseline0, reg = lm, y = 'scale(nbhsize)', xs = xs1) 
b <- reg.sum.bs(ls = cw_Baseline0, reg = lm, y = 'scale(outdegree)', xs = xs1)
a
b

aa <- reg.sum.bs(ls = cw_Baseline0, reg = lm, model = model3_sa) 
aa
aaa <- reg.sum.bs(ls = cw_Baseline0, reg = lm, model = model3_nb) 
aaa
## age and stage


####---- all lm ----
lapply(cw_Baseline0, function(x) {reg.sum.bs(ls = x, reg = lm, model = model3)
})
lapply(models, function(x) {reg.sum.bs(ls = l_EqualStage0, reg = lm, model = x)
})

####---- logistic ----
reg.sum.bs(ls = listUKclus, reg = glm, model = logit_model_uk, family = binomial(link = "logit"))

####---- stop ----