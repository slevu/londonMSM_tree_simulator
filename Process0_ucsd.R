####---- run model0.R ----####
## should run model with scenarios

#### rm(list=ls())

####---- lib ----
library(ape)
library(ggplot2)

if(TRUE){
  ####---- load sim ----
  ## Load newest Rdata
  l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
  load(l[length(l)])
  ls()
  simtree <- tree 
  rm(tree)

  ####---- load uk stuff ----
  t_uk <- read.tree(file = "../phylo-uk/data/ExaML_result.subUKogC_noDRM.finaltree.000")
  ## drop OG
  og <- c("Ref1", "Ref2", "Ref3", "Ref4", "Ref5", "Ref6", "HXB2")
  uktree <- drop.tip(t_uk, og ) 
 }

##-- Create edge list of (transformed) distances
if(FALSE){
  
  source("TReeToEdgeList.R") # returns path to RDS file
  
  ##- simtree
  system.time(
  path.simtree_el <- TreeToEdgeList(simtree, 
                       rate = 4.3e-3/365) # Berry et al. JVI 2007
  ) 
  
  ##- uktree
  system.time(
    path.uktree_el <- TreeToEdgeList(uktree, 
                        rate = 1) ) 
}

    
# Alternatively take TN93 distances
#  dukTN93 <- readRDS(file = "../phylo-uk/source/subUKogC_noDRM_151202_ucsdTN93.rds" )
  


####---- cluster UCSD ----
# source("TreeToClust.R")
source("EdgeListToClust.R")

## get list of quantiles, commands and list of
## dataframes of cluster members 
if(FALSE){
  system.time(
    simclus <- ucsd_hivclust(path.simtree_el, 
                             thr = c(0.015, 0.02, 0.05, 0.1))
  )
  
  system.time(
    ukclus <- ucsd_hivclust(path.uktree_el,
                            c(0.015, 0.02, 0.05, 0.1))
  )
  names(simclus)
  names(ukclus)
  # unlist(ukclus["warn"] , recursive = T, use.names = F)
  
  ###--- save 
  saveRDS(simclus, file = "data/simclus2.rds")
  saveRDS(ukclus, file = "data/ukclus2.rds")
 }

# rm(dsimtree, duktree)


####---- read ----
### only cluster members (based on given threshold)
simclus <- readRDS(file = "data/simclus2.rds")[["cl"]]
ukclus <- readRDS(file = "data/ukclus2.rds")[["cl"]]

### only cluster members (based on quantiles)
# simclus <- readRDS(file = "data/simclus.rds")[[3]]
# ukclus <- readRDS(file = "data/ukclus.rds")[[3]]

####---- quantiles ----
## read saved results of UCSD clustering
# readRDS(file = "data/simclus.rds")[[1]]
# readRDS(file = "data/ukclus.rds")[[1]]
####---- stop ----

####---- desc ----
##- Calculate size(=Freq) of each cluster across different threshold
simfreqClust <- lapply(simclus, 
                       function(x) as.data.frame(table(x$ClusterID), 
                       stringsAsFactors = FALSE))
ukfreqClust <- lapply(ukclus, 
                      function(x) as.data.frame(table(x$ClusterID), 
                      stringsAsFactors = FALSE))

##- number of different clusters by threshold
sapply(simfreqClust, function(x) dim(x)[1])
sapply(ukfreqClust, function(x) dim(x)[1])
##- cluster size
sapply(simfreqClust, function(x) summary(x$Freq))
sapply(ukfreqClust, function(x) summary(x$Freq))


####---- plot cluster size ----
# ## how many plots
a <- length(simfreqClust)
b <- length(ukfreqClust)

##- distr of cluster sizes: log(x) and Y untransformed
par(mfcol=c(2, max(a, b)))
for (i in 1:max(a, b)){
  h <- hist(log(ukfreqClust[[i]]$Freq), 
       main = paste("uk", names(ukfreqClust)[i]),
       xlab = "log(size)")
  hist(log(simfreqClust[[i]]$Freq), 
       main = paste("sim", names(simfreqClust)[i]),
       xlab = "log(size)")
}

####---- plot log-log ----
##- distr of cluster sizes: log(x) and log(y)
## how many plots
a <- length(simfreqClust)
b <- length(ukfreqClust)

par(mfcol=c(2, max(a, b)))
for (i in 1:max(a, b)){
  h <- hist(log(ukfreqClust[[i]]$Freq), plot = F)
  h$counts <- log1p(h$counts) # log(y)
  plot(h, ylab = "log(Freq)", 
       main = paste("uk", names(ukfreqClust)[i]),
       xlab = "log(size)")
  
  h <- hist(log(simfreqClust[[i]]$Freq), plot = F)
  h$counts <- log1p(h$counts) # log(y)
  plot(h, ylab = "log(Freq)",
          main = paste("sim", names(simfreqClust)[i]),
          xlab = "log(size)")
  
}

####---- QQ plot ----
par(mfcol=c(2, max(a, b)))
for (i in 1:max(a, b)){
  qqplot(ukfreqClust[[i]]$Freq, 
         simfreqClust[[i]]$Freq,
         main = paste("uk:", names(ukfreqClust)[i],
                      "sim:", names(simfreqClust)[i]),
         xlab = "uk", ylab = "sim")

  qqplot(log(ukfreqClust[[i]]$Freq), 
         log(simfreqClust[[i]]$Freq),
         main = paste("uk:", names(ukfreqClust)[i],
                      "sim:", names(simfreqClust)[i]),
         xlab = "log(uk)", ylab = "log(sim)")

}
dev.off()

####---- data ----
##- converting sample states in table of co-variates ?

if(FALSE){
# sampleTimes <- tree$sampleTimes
# sampleStates  <- tree$sampleStates 

demo <- matrix(NA, nrow =  length(tree$sampleTimes), ncol = 6)
for (i in 1:dim(tree$sampleStates)[1]){ # dim(ss)[1]
  deme <- names(which.max(tree$sampleStates[i,])) # name of column which has max value
  patient <- as.numeric(rownames(tree$sampleStates)[i])
  time <- tree$sampleTimes[i]
  age <- as.numeric( regmatches( deme, 
                regexec( "\\.age([0-9])", deme) )[[1]][2] )
  care <- as.numeric( regmatches( deme, 
                regexec( "care([0-9])", deme) )[[1]][2] )
  stage <- as.numeric( regmatches( deme, 
                regexec( "stage([0-9])", deme) )[[1]][2] )
  risk <- as.numeric( regmatches( deme, 
                regexec( "riskLevel([0-9])", deme) )[[1]][2] )
  demo[i,] <- cbind(patient, time, age, care, stage, risk)
} 
colnames(demo) <- cbind("patient", "time",
                        "age", "care", "stage", "risk")

saveRDS(as.data.frame(demo), file = "demo.rds")
}
####---- demo ----
demo <- readRDS("demo.rds")

##- date of diagnosis ?
date0 <- as.Date('1979-01-01')
demo$datediag <- date0 + demo$time
# min(demo$datediag) head(demo)
# max(demo$datediag)

###--- add cluster sizes
####... and outdegrees

### --- start function clust.stats --- ###
##- calculate both numclus and sizeclus for each seqindex into a LIST
##- with same variable names
##- For UCSD files, no ID if no cluster (size < 2) !

clust.stats <- function(clus = simclus, tree = simtree){
  
  ## get ALL tips names
  tip.names <- data.frame("id" = tree$tip.label , 
                          stringsAsFactors = F)
  ## get size of clusters
  freqClust <- lapply(clus, function(x){
    as.data.frame(table(x$ClusterID),
                  stringsAsFactors = FALSE)
  })
  
  ## empty list
  l <- list()
  
  ## loop over thresholds
  for (i in 1:length(clus) ) {
  ## merge cluster number (with NA)
   a <- merge(x = tip.names, y = clus[[i]],
          by.x = "id", by.y = "SequenceID",
          all.x = TRUE, sort = FALSE)
   
   ## merge cluster size (with NA)
   b <- merge(x = a, y = freqClust[[i]], 
              by.x = "ClusterID", by.y = "Var1", 
              all.x = TRUE, sort = FALSE)
   
  #- size 1 if not into a cluster
  b$Freq[is.na(b$Freq)] <- 1

  #- binary clustering variable
  b$binclus <- ifelse(b$Freq > 1 & !is.na(b$Freq), 1, 0)
  
  ##- pseudo ClusterID when size = 1
  ## starting from max(ClusterID+1)
  ## important for down-sampling
  .start <- max(b$ClusterID, na.rm = T)
  .n <- dim(b[is.na(b$ClusterID),] )[1]
   b$ClusterID[is.na(b$ClusterID)]  <- .start + 1:.n 

  #- colnames
  colnames(b)[which(colnames(b) =="Freq")] <- "size"
  l[[i]] <- b
  names(l)[i] <- names(freqClust[i])
  }
  return(l)
}
### --- end function clust.stats --- ###

l_sim <-  clust.stats(clus = simclus, tree = simtree)
l_uk <-  clust.stats(clus = ukclus, tree = uktree)
# str(l_sim)
# str(l_uk)

####---- proportion ----
##-proportion in or out clusters
sapply(l_sim, function(x) round(prop.table(table(x$binclus)),2))
sapply(l_uk, function(x) round(prop.table(table(x$binclus)),2))

##- cluster sizes (by individuals having such a size !!)
## probably meaningless
sapply(l_sim, function(x) summary(x$size))
sapply(l_uk, function(x) summary(x$size))

####---- merge ----
##- for simulated data
listclus <- lapply(l_sim, function(x) 
merge(x, demo, 
      by.x = "id", by.y = "patient", 
      all.x = T, sort = FALSE))

# head(listclus[[3]])
# table(listclus[[1]]$binclus, useNA = "ifany")
# table(listclus[[1]]$size, useNA = "ifany")

####--- sim naive regressions ---

####---- sim linear ----
###### just on low and high threshold
simli <- listclus[c(1:length(listclus))]

lm_model_factor <- "size ~ factor(stage) + factor(risk) + factor(age)"
lm_model_ordinal  <-  "scale(size) ~ scale(age) + scale(stage) + scale(time) + scale(risk)"
lm_model_ordinal_wo_time  <-  "scale(size) ~ scale(age) + scale(stage) + scale(risk)"

##- standard output
lapply(simli , function(x) summary(lm(lm_model_ordinal, data = x)))

##- return only parameters, p-values and R squared
lm.sum <- function(ls, lm_model){
  ## pvalue by threshold
  pvalue <- sapply(ls , function(x){
    coef(summary(lm(lm_model, data = x)))[,4]
  })
    ## coded by significance
    pvalue.code <- as.data.frame(
      apply(pvalue,2, function(x){
      cut(x,
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          labels=c('***','**', '*', '.', ''))
    }))
    row.names(pvalue.code) <- row.names(pvalue)
    ## parameter by threshold
  param <- signif(sapply(ls , function(x){
    coef(summary(lm(lm_model, data = x)))[,1]
  }), 2)
  ## R square
  r2 <- signif(sapply(ls , function(x){
    summary(lm(lm_model, data = x))$r.squared
  }), 3)
  return(list("parameter" = param, "pvalue" = pvalue.code, "r.squared" = r2))
}

lm.sum(ls = simli, lm_model = lm_model_ordinal)
lm.sum(ls = simli, lm_model = lm_model_factor)


##---- sim logistic ---- 
##- model: clus ~ age +  stage + time + risk
##- care = 1 for all at diagnosis
## ex. 
logit_model_std = "binclus ~ scale(age) + scale(stage) + scale(time) + scale(risk)"
lapply(simli , function(x) summary(glm(formula = logit_model_std,
                                   data = x,
                                   family = binomial(link = "logit"))
))

###--- Down-sample 1
###--- sort out dependency between indivduals from same cluster
### 1. downsample to make analysis of one cluster size
###  explained by median or mean of each co-variate.
###  2. plot the distribution of covariates by cluster size. 
###  With the intuition that smaller clusters are more explained
###  by covariates and larger ones are more random. 
###  Do it on real data and simulation

# For each cluster size, compute mean of all coavariates

####---- downsample ----
##- 1. down-sample: mean of each variable
down <- lapply(simli, function(x) aggregate(x[, 5:9], list("size" = x$size), mean))
# str(down) 
# head(down[[1]])
# head(simli[[1]])

##- linear regression
# lapply(down, function(x) summary(lm(lm_model_ordinal, data = x)))
lm.sum(ls = down, lm_model = lm_model_ordinal)

####---- lattice ----
##- 2. plots
library(lattice)
# trellis.par.set(canonical.theme(color = FALSE))
for(i in 1:length(simli)){
  print(histogram(~ stage|factor(size), 
            main = paste("distribution of stages by cluster sizes at threshold =", names(simli[i])),
            data = simli[[i]])
  )
}


##- plot correlation
##- df containing indepvar in a list
size.vs.covar <- function(l, depvar = "size",
                          indepvar = c("stage", "risk", "age", "time")){
  
  par(mfcol=c(length(indepvar), length(l) ), 
      mar = c(4,3,3,2)+0.1, oma = c(0, 0, 2, 0), bty = 'n') # b,l,t,r
  for(c in 1:length(l)){
    for(r in 1:length(indepvar)){
      plot( x = l[[c]][, indepvar[r]], y = l[[c]][, depvar], 
            ylab = '', xlab = '', font.main = 1,
            col="#00000050",
            main = paste(indepvar[r], names(down)[c]))
    }
  }
  mtext(paste(depvar, "(y) vs co-variates (x) "), outer = TRUE, cex = 1)
}

size.vs.covar(down)
dev.off()


## For ggplot, unlist ...

###--- Down-sample 2
### de-correlating: sample one individual by cluster. Repeat many times. Ensure much more power than downsampling with just one value per sample size.
head(listclus[[1]])
summary(listclus[[1]]$size) ## contains size = 1

###--- start function by threshold 

# df <- head(listclus[[1]],100)

downsample <- function(df, lm_model = lm_model, iter = 2){
  ##- sampling one id per cluster k times
  k <- iter
  ## loop
    ## empty list
    down_listclus <- vector("list", k)
    ##  k selection of one id from df
    ##   sampled by each ClusterID
    for (i in 1:k){
      down_listclus[[i]] <- df[df$id %in% 
                                 tapply(df$id,
                  df$ClusterID, 
                  function(x) sample(x, 1)),]
    }
  # head(down_listclus[[1]])
  # dim(down_listclus[[1]])
  
  ##- linear regression for each iteration
  fit <- lapply(down_listclus,
                function(x)
                  summary(lm(lm_model, data = x))) 
  
  ##- extract coefficient
    ## number of variables + intercept
    nvar <- dim(coef(fit[[1]]))[1]
    
    ## empty matrix of coefficients
    coef_lm <- matrix(NA, nvar * k, 4, 
                  dimnames = list(
                    rep(rownames( coef(fit[[1]]) ), k),
                    colnames(coef(fit[[1]])))) 
    ## loop
    for (i in 1:length(fit)){
      fill <- (i-1)*nvar + 1:nvar
      coef_lm[ fill,  ] <- coef(fit[[i]])
    }
  
  ##- number of p-value < 0.05
  sum <- tapply(coef_lm[,4], rownames(coef_lm), 
                function(x) sum(x < 0.05) / length(x))
  
  ##- mean of co-variates over iterations
  # head(down_listclus[[1]])
    ## change structure
    bind_df <- do.call(rbind, down_listclus)
    ## mean by ClusterID
    mean.sample <- aggregate(bind_df[, 5:9],
                           list("ClusterID" = bind_df$ClusterID, 
                                "size" = bind_df$size),                                                  mean)
  
  return(list(# down_listclus = down_listclus,
              #fit = fit, coef_lm = coef_lm,
              percent.signif =  sum, mean.sample = mean.sample))
}

###--- run down-sample over thresholds
## categorical variables
dd_cat <- sapply( listclus, function(x) {
  downsample(df = x,
             lm_model = lm_model_factor,
             iter = 100)
    }
  )
  ## percent of signif paramater
  t(do.call(rbind, dd_cat["percent.signif", ] ))

## ordinal variables
  dd_ord <- sapply( listclus, function(x) {
    downsample(df = x,
               lm_model = lm_model_ordinal,
               iter = 100)
  }
  )
  ## percent of signif paramater
  t(do.call(rbind, dd_ord["percent.signif", ] ))
  
  ## mean of co-variates by cluster
  # str(dd["mean.sample",])
  mean.down <- dd_ord["mean.sample",]

  ##- linear regression ordinal
  lapply(mean.down, function(x) 
    summary(lm(lm_model_ordinal, data = x)))
  lm.sum(ls = mean.down, lm_model = lm_model_ordinal)
  
  ## plot
  size.vs.covar(mean.down)
  dev.off()

##################################-----------# LLLLLAAAAAAAA
##################################
  ##- distribution of cluster size by stage, etc.
  
  boxplot(size ~ factor(stage), 
          main = 'size by stage', 
          data = simli[[4]])
  

## llllllllllllllllllaaaaaaaaaaa----------------------------> ??????
## make quantiles( .1, 1, 5 % of p-value)
## for each independant variable
## use test_lmTotable.R
## make plots of association
## loop at different thershold


###-------------------###
###--- for uk data ---###
###-------------------###

##---- multivariate ----
##- add demo outcome: stage of infection, treatment status, age group, CHIC or not
load("../phylo-uk/data/sub.RData")
rm(s)
##- selection of df covariates
y <- df[,c("seqindex","patientindex",  
           "agediag", "cd4", "vl", "onartflag",
           "ydiag", "agediag_cut", "cd4cut",
           "ydiag_cut", "CHICflag", "status")]

y$logvl <- log(y$vl)
y$sqrtcd4 <- sqrt(y$cd4)

#### get tips labels
# sim.names <- data.frame("id" = labels(dsimtree), 
#                         stringsAsFactors = F)
# saveRDS(sim.names, file = "sim.names.rds")
uk.names <- readRDS("uk.names.rds") 
##- in list 
l <- list()
for ( i in 1:length(ukclus) ) {
  ## merge cluster number (with NA)
  a <- merge(x = uk.names, y = ukclus[[i]],
             by.x = "id", by.y = "SequenceID",
             all.x = TRUE, sort = FALSE)
  
  ## merge cluster size (with NA)
  b <- merge(x = a, y = ukfreqClust[[i]], 
             by.x = "ClusterID", by.y = "Var1", 
             all.x = TRUE, sort = FALSE)
  
  #- size 1 if not into a cluster
  b$Freq[is.na(b$Freq)] <- 1
  
  #- binary clustering variable
  b$binclus <- ifelse(b$Freq > 1 & !is.na(b$Freq), 1, 0)
  
  #- colnames
  colnames(b)[which(colnames(b) =="Freq")] <- "size"
  l[[i]] <- b
  names(l)[i] <- names(ukfreqClust[i])
}
rm(a, b)
# str(l)

####---- proportion UK ----
##-proportion in or out clusters
sapply(l, function(x) round(prop.table(table(x$binclus)),2))
##- cluster sizes (by individuals having such a size !!)
sapply(l, function(x) summary(x$size))

####---- merge UK ----
listUKclus <- lapply(l, function(x) 
  merge(x, y, 
        by.x = "id", by.y = "seqindex", 
        all.x = T, sort = FALSE))

####--- naive regressions ---


####---- linear UK----
#### just on low and high threshold (but not too high !)
li <- listUKclus[ 1:(length(listUKclus)-1) ]
lm_model_std = "scale(size) ~ scale(agediag) + scale(sqrt(cd4)) +  scale(ydiag)"
lapply(li, function(x) summary(lm(lm_model_std, data = x)))

##---- logistic UK ---- 
##- model: clus ~ age +  stage + time + risk
##- care = 1 for all at diagnosis
## ex. 
logit_model_std = "binclus ~ scale(agediag) + scale(sqrt(cd4)) + scale(ydiag)"
lapply(li, function(x) summary(glm(formula = logit_model_std,
                                     data = x,
                                     family = binomial(link = "logit"))
                               ))


####---- downsample UK ----
##- 1. down-sample: MEDIAN of each variable (with na.rm = T)

# head(li[[1]][, c("agediag", "sqrtcd4", "ydiag")])
down_median <- lapply(li, function(x) 
  aggregate(x[, c("agediag", "sqrtcd4", "ydiag", "logvl")],
            list("size" = x$size), 
            function(x) median(x, na.rm = TRUE)))

down_mean <- lapply(li, function(x) 
  aggregate(x[, c("agediag", "sqrtcd4", "ydiag", "logvl")],
            list("size" = x$size), 
            function(x) mean(x, na.rm = TRUE)))
# str(down_mean) 
##- linear regression
lm_model_std = "scale(size) ~ scale(agediag) + scale(sqrtcd4) + scale(ydiag)"
lapply(down_median, function(x) summary(lm(lm_model_std, data = x)))
lapply(down_mean, function(x) summary(lm(lm_model_std, data = x)))

####---- lattice UK ----
##- 2. plots
library(lattice)
# trellis.par.set(canonical.theme(color = FALSE))
for(i in 1:length(li)){
  print(histogram(~ sqrtcd4|factor(size), 
                  main = paste("distribution of sqrt(cd4) by cluster sizes at threshold =", names(li[i])),
                  data = li[[i]])
  )
}

plot(x = li[[1]]$size, y = li[[1]]$sqrtcd4)

####---- surplus ----




####---- end ----

##-- check assortativity 
## first use EdgeList at each level
## then function AssortMix

##- false discovery rate: trnamission rates are equal according to cluster or SA ?
## on 100 trees, find significant association
## 
## study variance in cluster size, outdegrees changes with respect to the risk levels

##- UCSD cluster
##

########------ surplus


# ###--- bind table in one dataframe
# # list.files(tempdir())
# i <- 2
# chop <- function(patrn){
#   ## list outputs
#   l <- list.files(getwd(), pattern = patrn)
#   ## bind outputs
#   cl <- list()
#   for(i in 1:length(l)){
#     cl[[i]] <-  read.csv( l[i] )
#     # name threshold
#     names(cl)[i] <- substr(l[i], 
#                            regexpr(l[i], pattern = ".csv")[1] - 4,
#                            regexpr(l[i], pattern = ".csv")[1] - 1)
#   }
#   return(cl)
# }
# 
# simclus <- chop(patrn = "dsimtree")
# ukclus <- chop(patrn = "duktree")
# # str(simclus)
#  str(cl)
#  
# ####---- cluster UPGMA ----
# simhc <- hclust(dsimtree, method = "average") # UPGMA
# ukhc <- hclust(duktree, method = "average")
# 
# ##- cut based on height
# if (F){
# #- function of heights
# nheights <- 10 # number of threshold
# up <- round(mean(dsimtree), 1) # cut up to the mean
# # breaks <-  seq(1/nbreaks, 1-(1/nbreaks), by = 1/nbreaks)
# simclus <- cutree(simhc, h = seq(up / nheights, up, by = up / nheights) ) # h = breaks
# up <- round(mean(duktree),1) # cut up to the mean
# ukclus <- cutree(ukhc,  h = seq(up / nheights, up, by = up /nheights) )
# # rm(simx, ukx)
# }
# 
# #- cut as function of k groups
# kgroups <- c(500, 1000, 5000, 8000, 10000)
# simclus <- cutree(simhc, k = kgroups ) # h = breaks
# ukclus <- cutree(ukhc,  k = kgroups )
# # colnames(simclus) <- paste("k",colnames(simclus),sep='')
# # colnames(ukclus) <- paste("k",colnames(ukclus),sep='')
# 
# ###--- end UPGMA
