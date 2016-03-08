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

  ####---- load uk stuff ----
  t_uk <- read.tree(file = "../phylo-uk/data/ExaML_result.subUKogC_noDRM.finaltree.000")
  ## drop OG
  og <- c("Ref1", "Ref2", "Ref3", "Ref4", "Ref5", "Ref6", "HXB2")
  t_uk <- drop.tip(t_uk, og ) 
  t_uk
}

##-- Create edge list of transformed distances
source("TReeToEdgeList.R")
system.time(
 path.el <-  TreeToEdgeList(simtree)
)
# path.el <- "data/simtree_el.rds"
el <- readRDS(file = path.el)
head(el)




  # uk tree OR get TN93 distances ?
#   if (file.exists("data/uktree_dist.rds")){
#     duktree <- readRDS("data/uktree_dist.rds")
#   } else {
#     duktree <- as.dist(cophenetic.phylo(t_uk))
#     saveRDS(duktree, file = "data/uktree_dist.rds")
#   }
#  dukTN93 <- readRDS(file = "../phylo-uk/source/subUKogC_noDRM_151202_ucsdTN93.rds" )
  
  ## get and keep tips labels
#   if(FALSE){
#      sim.names <- data.frame("id" = labels(dsimtree), 
#                            stringsAsFactors = F)
#    saveRDS(sim.names, file = "sim.names.rds")
#    uk.names <- data.frame("id" = labels(duktree), 
#                            stringsAsFactors = F)
#    saveRDS(uk.names, file = "uk.names.rds")
#   }
  

#  hist(duktree, breaks = 50, xlab = "distance", ylab = "frequency", main = "UK MSM tree's distances", col = "grey")


####---- cluster UCSD ----
# source("TreeToClust.R")
source("EdgeListToClust.R")

## test for debug
## str(dsimtree)
## d <-  dsimtree 
# m <- as.matrix(dsimtree)
# d <- as.dist(m[1:200, 1:200])
# simclus <- ucsd_hivclust(d)
# str(d)
# head(cl)

## get list of quantiles, commands and list of
## dataframes of cluster members 
if(FALSE){
  system.time(
    simclus <- ucsd_hivclust(el)
  )
  
#   system.time(
#     ukclus <- ucsd_hivclust(duktree)
#   )
  
  ###--- save 
  saveRDS(simclus, file = "data/simclus.rds")
  saveRDS(ukclus, file = "data/ukclus.rds")
 }

# rm(dsimtree, duktree)


####---- read ----
### only cluster members
simclus <- readRDS(file = "data/simclus.rds")[[3]]
ukclus <- readRDS(file = "data/ukclus.rds")[[3]]

####---- quantiles ----
## read saved results of UCSD clustering
readRDS(file = "data/simclus.rds")[[1]]
readRDS(file = "data/ukclus.rds")[[1]]
####---- stop ----


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

###--- loop for
## For simulations,
##- function to calculate both numclus and sizeclus for each seqindex into a LIST
##- with same variable names
##- For UCSD files, no ID if no cluster (size < 2) !!!

## get tips labels
sim.names <- readRDS("sim.names.rds")
  ##- in list
  l <- list()
  for (i in 1:length(simclus)) {
  ## merge cluster number (with NA)
   a <- merge(x = sim.names, y = simclus[[i]],
          by.x = "id", by.y = "SequenceID",
          all.x = TRUE, sort = FALSE)
   
   ## merge cluster size (with NA)
   b <- merge(x = a, y = simfreqClust[[i]], 
              by.x = "ClusterID", by.y = "Var1", 
              all.x = TRUE, sort = FALSE)
   
  #- size 1 if not into a cluster
  b$Freq[is.na(b$Freq)] <- 1

  #- binary clustering variable
  b$binclus <- ifelse(b$Freq > 1 & !is.na(b$Freq), 1, 0)
  
  #- colnames
  colnames(b)[which(colnames(b) =="Freq")] <- "size"
  l[[i]] <- b
  names(l)[i] <- names(simfreqClust[i])
  }
  rm(a, b)
# str(l)

####---- proportion ----
##-proportion in or out clusters
sapply(l, function(x) round(prop.table(table(x$binclus)),2))
##- cluster sizes (by individuals having such a size !!)
sapply(l, function(x) summary(x$size))

####---- merge ----
listclus <- lapply(l, function(x) 
merge(x, demo, 
      by.x = "id", by.y = "patient", 
      all.x = T, sort = FALSE))

# head(listclus[[3]])
# table(listclus[[1]]$binclus, useNA = "ifany")

####--- sim naive regressions ---

####---- sim linear ----
###### just on low and high threshold
simli <- listclus[c(1:length(listclus))]

# lm_model_std = "size ~ age + stage + time + risk"
lm_model_w_time = "scale(size) ~ scale(age) + scale(stage) + scale(risk)"
lm_model_std = "scale(size) ~ scale(age) + scale(stage) + scale(time) + scale(risk)"
lapply(simli , function(x) summary(lm(lm_model_std, data = x)))

##- univariate
# summary(lm(
#   scale(size) ~ scale(stage),
#   data = simli[[3]]
# ))
# 
# summary(lm(
#   scale(size) ~ scale(risk),
#   data = simli[[1]]
# ))

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
###--- sort dependency between indivduals from same cluster
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

##- linear regression
lapply(down, function(x) summary(lm(lm_model_std, data = x)))

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

###--- Down-sample 2
### de-correlating: sample one individual by cluster. Repeat many times. Ensure much more power than downsampling with just one value per sample size.
head(listclus[[1]])

## petit
#df <- head(listclus[[1]], 10)
## grand
lm_model <- lm_model_w_time # lm_model_std
df <- listclus[[1]]
## sampling one id per cluster k times
k <- 100
## empty list
down_listclus <- vector("list", k)
## loop: k selection of id from df sampled by ClusterID
for (i in 1:k){
  down_listclus[[i]] <- df[df$id %in% tapply(df$id,
                  df$ClusterID, 
                  function(x) sample(x, 1)),]
                  
}
dim(down_listclus[[1]])[1]

##- linear regression
fit <- lapply(down_listclus, 
       function(x) 
         summary(lm(lm_model, data = x))) 
## extract coefficient
length(fit)

## test i <- 1
# coef_lm <- vector("list", k)
# for (i in 1:length(fit)){
#   a <- (capture.output(fit[[i]]))
# #   grep("Intercept", a)
# #   grep("---", a)
#   b <- a[(grep("Intercept", a)):(grep("---", a)-1)]
#   
#   coef_lm[[i]] <- read.table(textConnection( b ), fill = TRUE)
#   names(coef_lm[[i]]) <- c("Var", colnames(coef(fit[[1]]) ), "Signif")
# }
# coef_lm

## 2nd try

## try to allocate the number of variables + intercept
nvar <-  ifelse( gregexpr("\\+", lm_model)[[1]][1] == -1, 2,
        length( gregexpr("\\+", lm_model)[[1]] ) + 2
        ) 

## matrix of coefficients
coef_lm <- matrix(NA, nvar * k, 4, 
                  dimnames = list(
                    rep(rownames( coef(fit[[1]]) ), k),
                    colnames(coef(fit[[1]])))) 
# i <- 1
for (i in 1:length(fit)){
  fill <- (i-1)*nvar + 1:nvar
  coef_lm[ fill,  ] <- coef(fit[[i]])
}
## number of p-value < 0.05
# coef_lm[rownames(coef_lm) == "scale(stage)", 4]
tapply(coef_lm[,4], rownames(coef_lm), function(x) sum(x < 0.05))
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
