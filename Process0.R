####---- run model0.R ----####
## should run model with scenarios

#### rm(list=ls())

####---- lib ----
library(ape)
library(ggplot2)

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

####---- get distances ----
####  cluster size to real data. Need to have same number of clusters ?

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
head(dsimtree)
head(duktree)

## normalize
simx <- dsimtree / (max(dsimtree) - min(dsimtree))
ukx <- duktree / (max(duktree) - min(duktree))
# rm(simtree, uktree)

##- histogram distances
# summary(x)
hist(simx, breaks = 50, xlab = "distance", ylab = "frequency", main = "", col = "grey")
hist(ukx, breaks = 50, xlab = "distance", ylab = "frequency", main = "", col = "grey")

####---- cluster UPGMA ----
simhc <- hclust(simx, method = "average") # UPGMA
ukhc <- hclust(ukx, method = "average")

##- cut based on height
if (F){
#- function of heights
nheights <- 10 # number of threshold
up <- round(mean(simx), 1) # cut up to the mean
# breaks <-  seq(1/nbreaks, 1-(1/nbreaks), by = 1/nbreaks)
simclus <- cutree(simhc, h = seq(up / nheights, up, by = up / nheights) ) # h = breaks
up <- round(mean(ukx),1) # cut up to the mean
ukclus <- cutree(ukhc,  h = seq(up / nheights, up, by = up /nheights) )
# rm(simx, ukx)
}

#- cut as function of k groups
kgroups <- c(500, 1000, 5000, 8000, 10000)
simclus <- cutree(simhc, k = kgroups ) # h = breaks
ukclus <- cutree(ukhc,  k = kgroups )
# colnames(simclus) <- paste("k",colnames(simclus),sep='')
# colnames(ukclus) <- paste("k",colnames(ukclus),sep='')
head(simclus)
head(ukclus)

##- Calculate size(=Freq) of each cluster across different threshold
simfreqClust <- apply(simclus, 2, function(x) as.data.frame(table(x))) # list
ukfreqClust <- apply(ukclus, 2, function(x) as.data.frame(table(x)))
# str(simfreqClust)
# head(simfreqClust[[1]])

##- number of different clusters by threshold # if number varies !
# sapply(simfreqClust, function(x) dim(x)[1])
# sapply(ukfreqClust, function(x) dim(x)[1])

##- cluster size
#- sim
sapply(simfreqClust, function(x) summary(x$Freq))
#- uk
sapply(ukfreqClust, function(x) summary(x$Freq))
##- percentiles
# sapply(freqClust, function(x) round(quantile(x$Freq, 
# probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 0.99, 1))))
# 

####---- plot cluster size ----
##- distr of cluster sizes: log(x) and Y untransformed
par(mfcol=c(2, length(kgroups)))
for (i in 1:length(kgroups)){
  h <- hist(log(ukfreqClust[[i]]$Freq), 
       main = paste("uk", names(ukfreqClust)[i]),
       xlab = "log(size)")
  hist(log(simfreqClust[[i]]$Freq), 
       main = paste("sim", names(simfreqClust)[i]),
       xlab = "log(size)")
}

##- distr of cluster sizes: log(x) and log(y)
par(mfcol=c(2, length(kgroups)))
for (i in 1:length(kgroups)){
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
par(mfcol=c(2, length(kgroups)))
for (i in 1:length(kgroups)){
  qqplot(ukfreqClust[[i]]$Freq, 
         simfreqClust[[i]]$Freq,
         main = names(ukfreqClust)[i],
         xlab = "uk", ylab = "sim")

  qqplot(log(ukfreqClust[[i]]$Freq), 
         log(simfreqClust[[i]]$Freq),
         main = names(ukfreqClust)[i],
         xlab = "log(uk)", ylab = "log(sim)")

}

####---- ggplot cluster size ----
##- ggplot: distributions of cluster sizes for different number of clusters
##- dataframe: 
df <- data.frame()
for (i in 2:length(kgroups)){
  df <- rbind(df,
              rbind( cbind(gr = "sim", 
                           k = names(simfreqClust)[i], 
                           simfreqClust[[i]]) ,
                     cbind(gr = "uk", 
                           k = names(ukfreqClust)[i],
                           ukfreqClust[[i]]) )
              )
}
head(df)
str(df)


##- UPGMA comparison of cluster sizes
g <- ggplot(data = df, aes(Freq))
g +  geom_histogram(origin = 0, binwidth = .5) + 
  scale_x_continuous(breaks=c(0,1,3,10,30,100,300,1000,5000), 
                     trans = "log1p",
                     name = "Cluster size") +
  scale_y_continuous(breaks=c(0, 30, 100, 300, 1000),
                     trans = "log1p",
                     name = "Number of clusters") +
  facet_grid(k  ~ gr, 
             scales = "free",
             space = "free") +
  theme_minimal()
 # theme_bw()


####---- data ----
##- converting sample states in table of co-variates ?
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
str(demo)

##- date of diagnosis ?
date0 <- as.Date('1979-01-01')
demo$datediag <- date0 + demo$time
min(demo$datediag)
max(demo$datediag)

####---- add cluster sizes ----
####... and outdegrees

##---- loop ----
##- function to calculate both numclus and sizeclus for each seqindex into a LIST
##- with same variable names
  ##- in list
  l <- list()
  for (i in 1:length(kgroups)) {
    
  #- cluster number
  numclus <- as.data.frame(simclus[, i])
  numclus <- cbind(rownames(numclus), numclus)
  colnames(numclus) <- c("id", "num")
  row.names(numclus) <- NULL
  # head(numclus)
  
  #- size of cluster
  a <- merge(x = numclus, y = simfreqClust[[i]], 
             by.x = "num", by.y = "x", 
             all.x = TRUE, sort = FALSE)
  #- binary clustering variable
  a$Clus <- ifelse(a$Freq > 1, 1, 0)
  #- colnames
  colnames(a)[which(colnames(a) =="Freq")] <- "size"
  colnames(a)[which(colnames(a) =="Clus")] <- "clus"
  l[[i]] <- a
  names(l)[i] <- names(simfreqClust[i])
  }

  rm(a, numclus)
# str(l)

##-proportion in or out clusters
sapply(l, function(x) round(prop.table(table(x$clus)),2))
##- cluster sizes
sapply(l, function(x) summary(x$size))

##---- merge ----
listclus <- lapply(l, function(x) 
merge(x, demo, 
      by.x = "id", by.y = "patient", 
      all.x = T, sort = FALSE))

# head(listclus[[4]])
# table(listclus[[4]]$clus)

##---- logistic ---- 
##- model: clus ~ age +  stage + time + risk
##- care = 1 for all at diagnosis
## ex. 
logit_model = "clus ~ age + stage + time + risk"
logit_model_std = "clus ~ scale(age) + scale(stage) + scale(time) + scale(risk)"
?glm
lapply(listclus, function(x) summary(glm(formula = logit_model_std,
                                 data = x,
                                 family = binomial(link = "logit"))))

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
##
y <- rt(200, df = 5)
qqnorm(y); qqline(y, col = 2)
qqplot(y, rt(300, df = 5))

qqnorm(precip, ylab = "Precipitation [in/yr] for 70 US cities")