####---- run model0.R ----####
## should run model with scenarios

#### rm(list=ls())

####---- include ----
detail_knitr <- FALSE
source("functions.R")

####---- lib ----
library(ape)

if(TRUE){
  ####---- load sim ----
  ## Load newest Rdata
  l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
  l
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

####---- Create edge list ----
## of (transformed) distances
if(FALSE){
  
  # function TReeToEdgeList returns path to RDS file
  
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

####---- desc2 ----
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

####---- clust.stats ----
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
# saveRDS(listclus, file = "data/listclus_sim.rds")

###########################
### --- Regressions --- ###
###########################

####---- models ----
##- linear
lm_model_factor <- "size ~ factor(stage) + factor(risk) + factor(age)"
lm_model_ordinal  <-  "scale(size) ~ scale(age) + scale(stage) + scale(time) + scale(risk)"
lm_model_ordinal_wo_time  <-  "scale(size) ~ scale(age) + scale(stage) + scale(risk)"

##- logistic
logit_model_ord  <-  "binclus ~ scale(age) + scale(stage) + scale(time) + scale(risk)"
logit_model_fact <-  "binclus ~ factor(stage) + factor(risk) + factor(age)"

##- standard output
lapply(listclus , function(x) summary(lm(lm_model_ordinal, data = x)))

####--- sim naive regressions ---
####---- summary naive ----
reg.sum(ls = listclus, reg = lm, model = lm_model_ordinal)
reg.sum(ls = listclus, reg = lm, model = lm_model_factor)

##- logistic
##- model: clus ~ age +  stage + time + risk
##- care = 1 for all at diagnosis

# logfit <- lapply(simli , function(x){
#   summary(glm(formula = logit_model_std, 
#               data = x, family = binomial(link = "logit")))
#   }) 
####---- sim logistic ---- 
reg.sum(ls = listclus, reg = glm, model = logit_model_ord, family = binomial(link = "logit"))

reg.sum(ls = listclus, reg = glm, model = logit_model_fact, family = binomial(link = "logit"))
####---- stop ----

####---- downsample2 ----
### De-correlating: sample one individual by cluster. Repeat many times. Ensure much more power than downsampling with just one value per sample size.
head(listclus[[1]])
summary(listclus[[1]]$size) ## contains size = 1

###--- start function by threshold 
# df <- listclus[[1]]

downsample <- function(df, lm_model = lm_model, var = c("time", "age", "care", "stage", "risk"), iter = 2){
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
    ## mean by ClusterID of set of variables var
    mean.sample <- aggregate(bind_df[, var],
                           list("ClusterID" = bind_df$ClusterID, 
                                "size" = bind_df$size),                                                  mean)
  
  return(list(# down_listclus = down_listclus,
              #fit = fit, coef_lm = coef_lm,
              percent.signif =  sum, mean.sample = mean.sample))
}

####---- run down-sample 1 ----
### over thresholds

## ordinal variables
  dd_ord <- sapply( listclus, function(x) {
    downsample(df = x,
               lm_model = lm_model_ordinal,
               var = c("time", "age", "care", "stage", "risk"),
               iter = 100)
  }
  )
  ## percent of signif paramater
  t(do.call(rbind, dd_ord["percent.signif", ] ))
 
## categorical variables
  dd_cat <- sapply( listclus, function(x) {
    downsample(df = x,
               lm_model = lm_model_factor,
               var = c("time", "age", "care", "stage", "risk"),
               iter = 100)
  }
  )
  ## percent of signif paramater
  t(do.call(rbind, dd_cat["percent.signif", ] ))
  
####---- run down-sample 2 ----
  
  ## mean of co-variates by cluster
  # str(dd["mean.sample",])
  mean.down <- dd_ord["mean.sample",]

  ##- linear regression ordinal
#    lapply(mean.down, function(x) {
#    summary(lm(lm_model_ordinal, data = x))})
  reg.sum(ls = mean.down, reg = lm, model = lm_model_ordinal)
  
  ## plot
  size.vs.covar(mean.down)
  # dev.off()


####---- STOP ----######

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

####---- merge UK ----
listUKclus <- lapply(l_uk, function(x) 
  merge(x, y, 
        by.x = "id", by.y = "seqindex", 
        all.x = T, sort = FALSE))
# head(listUKclus[[1]])

####--- naive regressions ---


####---- linear UK----
#### just on low and high threshold (but not too high !)
# li <- listUKclus[ 1:(length(listUKclus)-1) ]
lm_model_uk = "scale(size) ~ scale(agediag) + scale(sqrt(cd4)) +  scale(ydiag)"
# lapply(li, function(x) summary(lm(lm_model_std, data = x)))
reg.sum(ls = listUKclus, reg = lm, model = lm_model_uk)

##---- logistic UK ---- 
##- model: clus ~ age +  stage + time + risk
##- care = 1 for all at diagnosis
## ex. 
logit_model_uk = "binclus ~ scale(agediag) + scale(sqrt(cd4)) + scale(ydiag)"

reg.sum(ls = listUKclus, reg = glm, model = logit_model_uk, family = binomial(link = "logit"))


####---- downsample UK ----
### de-correlating: sample one individual by cluster. Repeat many times. Ensure much more power than downsampling with just one value per sample size.
head(listUKclus[[1]])
summary(listUKclus[[1]]$size) ## contains size = 1

####---- run down-sample UK 1 ----
### over thresholds

## ordinal variables
dd_ord_uk <- sapply( listUKclus, function(x) {
  downsample(df = x,
             lm_model = lm_model_uk,
             var = c("agediag", "cd4", "ydiag"),
             iter = 100)
}
)
## percent of signif paramater
t(do.call(rbind, dd_ord_uk["percent.signif", ] ))

####---- run down-sample UK 2 ----

## mean of co-variates by cluster
# str(dd["mean.sample",])
mean.down_uk <- dd_ord_uk["mean.sample",]
# head(mean.down_uk[[2]])
##- linear regression ordinal
#    lapply(mean.down, function(x) {
#    summary(lm(lm_model_ordinal, data = x))})
reg.sum(ls = mean.down_uk, reg = lm, model = lm_model_uk)

## plot
size.vs.covar(l = mean.down_uk, depvar = "size",
indepvar = c("agediag", "cd4", "ydiag"))
# dev.off()

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
