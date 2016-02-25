####---- run model0.R ----####
## should run model with scenarios

#### rm(list=ls())

####---- lib ----
# library(ape)
# library(ggplot2)

if(FALSE){
  ####---- load sim ----
  ## Load newest Rdata
  l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
  load(l[length(l)])
  ls()
  tree 
  
  ####---- load uk stuff ----
  t_uk <- read.tree(file = "../phylo-uk/data/ExaML_result.subUKogC_noDRM.finaltree.000")
  ## drop OG
  og <- c("Ref1", "Ref2", "Ref3", "Ref4", "Ref5", "Ref6", "HXB2")
  t_uk <- drop.tip(t_uk, og ) 
  t_uk
}

####---- get distances ----
####  cluster size to real data. Need to have same number of clusters ?
if(FALSE){
  ## get distances
  ##- matrix first into distances
  
  #- sim tree
  if (file.exists("data/simtree_dist.rds")){
    dsimtree <- readRDS("data/simtree_dist.rds")
  } else {
    dsimtree <- as.dist(cophenetic.phylo(tree))
    saveRDS(dsimtree, file = "data/simtree_dist.rds")
  }
  # uk tree
  if (file.exists("data/uktree_dist.rds")){
    duktree <- readRDS("data/uktree_dist.rds")
  } else {
    duktree <- as.dist(cophenetic.phylo(t_uk))
    saveRDS(duktree, file = "data/uktree_dist.rds")
  }
  
  ## get and keep tips labels
   sim.names <- data.frame("id" = labels(dsimtree), 
                           stringsAsFactors = F)
   saveRDS(sim.names, file = "sim.names.rds")
   uk.names <- data.frame("id" = labels(duktree), 
                           stringsAsFactors = F)
   saveRDS(uk.names, file = "uk.names.rds")
  
  # head(dsimtree)
  # head(duktree)
  
  ## normalize
  dsimtree <- dsimtree / max(dsimtree)
  duktree <- duktree / max(duktree)
  
  ##- histogram distances
  # summary(x)
  hist(dsimtree, breaks = 50, xlab = "distance", ylab = "frequency", main = "Simulated tree's distances", col = "grey")
  hist(duktree, breaks = 50, xlab = "distance", ylab = "frequency", main = "UK MSM tree's distances", col = "grey")
}
####---- cluster UCSD ----
source("TreeToClust.R")

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
    simclus <- ucsd_hivclust(dsimtree)
  )
  
  system.time(
    ukclus <- ucsd_hivclust(duktree)
  )
  
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
demes <- as.vector(read.csv(file = "demes.csv")$x)
sampleTimes <- scan( file = 'sampleTimes' )
ss  <- matrix( scan( file = 'sampleStates' ) , 
               byrow = TRUE, 
               ncol = length(demes))
colnames(ss) <- demes
dim(ss)
max(ss[,121]) # nothing on source

demo <- data.frame()
for (i in 1:dim(ss)[1]){ # dim(ss)[1]
  deme <- names(which(ss[i,] == 1)) # name of column which has value 1
  patient <- i
  time <- sampleTimes[i]
  age <- as.numeric( regmatches( deme, 
                regexec( "\\.age([0-9])", deme) )[[1]][2] )
  care <- as.numeric( regmatches( deme, 
                regexec( "care([0-9])", deme) )[[1]][2] )
  stage <- as.numeric( regmatches( deme, 
                regexec( "stage([0-9])", deme) )[[1]][2] )
  risk <- as.numeric( regmatches( deme, 
                regexec( "riskLevel([0-9])", deme) )[[1]][2] )
  demo <- rbind(demo, cbind(
    patient, time, age, care, stage, risk))
} 
# str(demo)
saveRDS(demo, file = "demo.rds")
}
####---- demo ----
demo <- readRDS("demo.rds")

##- date of diagnosis ?
date0 <- as.Date('1979-01-01')
demo$datediag <- date0 + demo$time
# min(demo$datediag)
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
print("coucou")
# head(listclus[[3]])
# table(listclus[[1]]$binclus, useNA = "ifany")

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
## just on low and high threshold
l <- listclus[c(1,length(listclus))]
down <- lapply(l, function(x) aggregate(x[, 5:9], list("size" = x$size), mean))
# str(down)
# 
##- linear regression
lm_model_std = "scale(size) ~ scale(age) + scale(stage) + scale(time) + scale(risk)"
lapply(down, function(x) summary(lm(lm_model_std, data = x)))

####---- lattice ----
##- 2. plots
library(lattice)
# trellis.par.set(canonical.theme(color = FALSE))
for(i in 1:length(l)){
  print(histogram(~ stage|factor(size), 
            main = paste("distribution of stages by cluster sizes at threshold =", names(l[i])),
            data = l[[i]])
  )
}


###--- for uk data ---


##---- multivariate ----
##- add demo outcome: stage of infection, treatment status, age group, CHIC or not
load("../phylo-uk/data/sub.RData")
##- selection of df covariates
y <- df[,c("seqindex","patientindex",  
           "agediag", "cd4", "vl", "onartflag",
           "ydiag", "agediag_cut", "cd4cut",
           "ydiag_cut", "CHICflag", "status")]

y$logvl <- log(y$vl)
y$logcd4 <- log(y$cd4)

#### get tips labels
# sim.names <- data.frame("id" = labels(dsimtree), 
#                         stringsAsFactors = F)
# saveRDS(sim.names, file = "sim.names.rds")
uk.names <- readRDS("uk.names.rds") 
##- in list
l <- list()
for (i in 1:length(ukclus)) {
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

####---- downsample UK ----
##- 1. down-sample: median of each variable (with na.rm = T)
## just on low and high threshold (but not too high !)
l <- listUKclus[c(1,length(listUKclus)-1)]
# head(l[[1]][, c("agediag", "logcd4", "ydiag")])
down <- lapply(l, function(x) 
  aggregate(x[, c("agediag", "logcd4", "ydiag", "logvl")],
            list("size" = x$size), 
            function(x) median(x, na.rm = TRUE)))
# str(down) 
# str(l)
# 
##- linear regression
lm_model_std = "scale(size) ~ scale(agediag) + scale(logcd4) + scale(logvl) + scale(ydiag)"
lapply(down, function(x) summary(lm(lm_model_std, data = x)))

####---- lattice UK ----
##- 2. plots
library(lattice)
# trellis.par.set(canonical.theme(color = FALSE))
for(i in 1:length(l)){
  print(histogram(~ logcd4|factor(size), 
                  main = paste("distribution of log(cd4) by cluster sizes at threshold =", names(l[i])),
                  data = l[[i]])
  )
}

plot(x = l[[1]]$size, y = l[[1]]$logcd4)

####---- surplus ----

##---- logistic ---- 
##- model: clus ~ age +  stage + time + risk
##- care = 1 for all at diagnosis
## ex. 
logit_model_std = "binclus ~ scale(agediag) + scale(logcd4) + scale(ydiag)"
lapply(l[1], function(x) summary(glm(formula = logit_model_std,
                                 data = l[[1]],
                                 family = binomial(link = "logit"))) 
       )

###################### missing values !!!!!!!! ###############

####---- surplus ----
# logistic <- function(x, m = logit_model){
#   fit <- glm(m , data = x, 
#                    family = binomial(link = "logit"))
#   co <- coef(summary(fit))
#   ## odds ratios and 95% CI
#   # or <- exp(cbind(OR = coef(fit), confint(fit)))
#   # return(list(co, or))
#   return(cbind(co[,c(1,4)]))
# }
# 
# ##- test 1 level
# c <- listclus[[1]]
# logistic(x = c, m = logit_model_std)
# 
# ##- all levels
# lapply(listclus, function(x) logistic(x, m = logit_model_std))

##---- linear ----
lm_model = "size ~ age + stage + time + risk"
lm_model_std = "scale(size) ~ scale(age) + scale(stage) + scale(time) + scale(risk)"

lapply(listclus, function(x) summary(lm(lm_model_std, data = x)))


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