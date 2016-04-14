###--- regressions of UK bootstrap ML trees, cluster and outdegree ---###
library(ape)
source("functions.R")

####---- include2 ----
detail_knitr <- TRUE


####---- stop ----
if(FALSE){
####---- get sampleTimes ----
##- list of lsd trees
## filename pattern from LSD changes with LSD version !
# list.lsd.trees <- list.files(path = "data/LSD", pattern = "result.date", full.names = TRUE) # macbook
list.lsd.trees <- list.files(path = "data/LSD", pattern = "result_newick_date", full.names = TRUE)
head(list.lsd.trees)

##- read first LSD tree to name 
## sampling times and CD4s with tip.labels
t <- read.tree( list.lsd.trees[2] )
STFN <- "data/LSD/t000.dates"

##- sampling times
dates <-( read.table(STFN, skip=1, 
                     colClasses=c('character', 'numeric') ) )
##- named vector
sampleTimes <- setNames( dates[,2], dates[,1] )[t$tip.label] 
head(sampleTimes)
rm(list.lsd.trees, STFN, t, dates)
}

####---- reprise ----
###--- Do these things to outdegree and cluster sizes :
## restrict or not to cohort of sampling 
## add explanatory variables
## categorize age and cd4
## run model and tests

####---- list of W ----
## list of infector prob files
list.W0 <- list.files("data/phydynR", pattern = 'mh20', full.names = TRUE)
## order
list.W0 <- list.W0[order(nchar(list.W0), list.W0)]

####---- analyses ----
##- calculating outdegree
##- but first restrict to cohort 
####---- time restriction ----
thr_year <- Inf

####---- list of outdegrees ----
if(FALSE){
#### for m bootstrap
### function: input filename of W
outdegree <- function(w.fn, t = thr_year){
  W <- readRDS(w.fn)
  
  ## restrict to cohort sampled within thr_year years
  cohort <- names(sampleTimes[sampleTimes > (max(sampleTimes) - t ) ] )
  i <- which( (W$donor %in% cohort) & (W$recip %in% cohort ) )
  WW <- list( donor = W$donor[i] , recip = W$recip[i], infectorProbability = W$infectorProbability[i] )
  
  ## calculate outdegrees
  out <- aggregate(
    x = list(outdegree = WW$infectorProbability),
    by = list(patient = WW$donor), 
    FUN = function(x) sum(x, na.rm = T) )
  return(out)
} 

list.outdegree <- lapply(list.W0, outdegree)
}
####---- stop ----
#  ## nb patients
#  summary(sapply(list.outdegree, function(x) dim(x)[1] ))
#  ## head
#  head(list.outdegree[[1]])
#  ## summary OD
#  summary(list.outdegree[[1]]$outdegree)

###--- clusters ---###
## get cluster size for one threshold (say 2%) in parallel with outdegree before doing regressions
####---- get cluster list ----
l_bs_uk <- readRDS( file = "data/listUK_ucsd_clus.rds")

####---- prune ----

##- function to prune cluster according to a threshold of sampling time to control for cohort effect
##- recalculate size and cluster membership
## or use ydiag ? which is different (median diff # 2.5 years)
## depends on clustering algorithm ?

prune.clus <- function(a, t = thr_year){
  ## subset df by sampling times
  cohort <- names(sampleTimes[sampleTimes > (max(sampleTimes) - t ) ] )
  aa <- a[a[,"id"] %in% cohort,]
  if(identical(aa,a)){ 
    print('do nothing')
    } else {
  ## for each clusterID, re-calculate size and binclus membership
  for (i in unique(aa[,"ClusterID"])){
    aa[ aa[,"ClusterID"] == i, "size" ] <- nrow(aa[ aa[,"ClusterID"] == i, ])
  }
  aa[,"binclus"] <- ifelse(aa[,"size"] < 2, 0, 1 )
    }
  return(aa)
}

list.clus.pruned <- lapply(l_bs_uk, function(x){
  lapply(x, prune.clus)
})

####---- demo variables ----
##- add individual explanatory variates
##- selection of df covariates
load("../phylo-uk/data/sub.RData")
y <- df[,c("seqindex","patientindex", 
           "dob_y", "agediag", "cd4", 
           "ydiag", "CHICflag", "ethnicityid")]

y$ethn.bin <- ifelse(y$ethnicityid == "White", "white", "not white")
y$CHICflag <- ifelse(y$CHICflag == "Yes", 1, 0)
y$ethnicityid <- NULL
y <- unfactorDataFrame(y)

## categorize continuous variables
y$agecl <- sapply( y[ , "agediag"] , age2quantile )
y$cd4cl <- sapply( y[ , "cd4"] , cd4toStage )
head(y)
rm(s, df)

####---- add demo ----
#### variables
if(FALSE){
cluster <- lapply(list.clus.pruned, function(u){
  lapply(u, function(x) {
    merge(x, y, 
        by.x = "id", by.y = "seqindex", 
        all.x = T, sort = FALSE)
})
  })

od <- lapply(list.outdegree, function(x){
  merge(x, y, 
        by.x = "patient", by.y = "seqindex", 
        all.x = T, sort = FALSE)
})
}
####---- load cluster and od ----
###- save and read
# saveRDS(cluster, file = "data/list_cluster_uk_bs_thr_demo.rds")
# saveRDS(od, file = "data/list_outdegree_uk_bs_demo.rds")
cluster <- readRDS(file = "data/list_cluster_uk_bs_thr_demo.rds")
od <- readRDS(file = "data/list_outdegree_uk_bs_demo.rds")

list.total <- c("SA" = list(od), "Cluster" = cluster)
names(list.total)

####---- stop ----
##- difference between sampling time and time of diagnosis
# both <- 
#   merge(a, 
#         data.frame("id" = names(sampleTimes), sampleTimes, row.names = NULL ), 
#         by = "id" )
# diff <- df$dateres - df$datediag
# summary(as.numeric(diff)/365)
# hist(both$sampleTimes - both$ydiag)
# summary(both$sampleTimes - both$ydiag)

## plot on first bootstrap
## 
# boxplot(sqrt(od[[1]]$outdegree) ~ od[[1]]$agecl)
# boxplot(sqrt(od[[1]]$outdegree) ~ od[[1]]$cd4cl)
# boxplot(sqrt(od[[1]]$outdegree) ~ od[[1]]$CHICflag)
# 
# boxplot(sqrt(cluster[[1]]$size) ~ cluster[[1]]$agecl)
# boxplot(sqrt(cluster[[1]]$size) ~ cluster[[1]]$cd4cl)
# boxplot(sqrt(cluster[[1]]$size) ~ cluster[[1]]$CHICflag)

############################################
###########-------- la ------############### 
############################################    

## make parallel results cluster - SA for R2, p-value, and effects

##- downsample ?

####---- fn ----
source("test_fn_compare.reg.sum.bs.R")
compare.reg.bs

####---- models ----
model1 <- "y ~ factor(agecl)"
model1c <- "y ~ scale(agediag)"
model1i <- "y ~ factor(agecl) + factor(cd4cl)"
model2 <- "y ~ factor(cd4cl)"
model3 <- "y ~ factor(ethn.bin)"
model4 <- "y ~ factor(CHICflag)"
model5 <- "y ~ factor(agecl) + factor(cd4cl) + factor(agecl)*factor(cd4cl) "
model6 <- "y ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn.bin) + factor(CHICflag)"

####---- tests ----

test1c <- compare.reg.bs(ls = list.total, reg = lm, model = model1c, alpha = 0.05)
test2 <- compare.reg.bs(ls = list.total, reg = lm, model = model2, alpha = 0.05)

####---- model age ----
test <- compare.reg.bs(ls = list.total, reg = lm, model = model1, alpha = 0.05)
test

####---- model age cd4 ----
test2 <- compare.reg.bs(ls = list.total, reg = lm, model = model1i, alpha = 0.05)
test2


####---- model continuous ----
test6 <- compare.reg.bs(ls = list.total, reg = lm, model = model6, alpha = 0.05)
test6

####---- model factor ----
test5 <- compare.reg.bs(ls = list.total, reg = lm, model = model5, alpha = 0.05)
test5


####---- end ----
model5 <- "y ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn.bin) + factor(CHICflag)"
model6 <- "y ~ factor(agecl) + factor(cd4cl) + factor(CHICflag) + factor(agecl)*factor(cd4cl) + factor(CHICflag)*factor(cd4cl)"
models <- as.list(paste0("model", 1:5))

####---- all lm ----
lapply(models, function(x) {compare.reg.bs(
  ls = list.total, reg = lm, model = x)
})
##- univariate
# lapply(models, function(x) {reg.sum.bs(ls = listUKclus, reg = lm, model = x)
# })
compare.reg <- function()
  
  
  y <- "scale(outdegree)"
full.model = paste(y, model6)
reg.sum.bs(ls = list(od), reg = lm, model = full.model)

y = "scale(size)"
full.model = paste(y, model6)
reg.sum.bs(ls = list(cluster), reg = lm, model = full.model) 

##- examples
# summary(lm("size ~ scale(CHICflag) + scale(sqrt(cd4)) ", cluster[[1]][[1]]))
# summary(
#   glm("binclus ~ factor(CHICflag) ",
#       family = binomial(link = "logit"), cluster[[3]])
#   )

###-- regressions
###- test interaction for CHICflag*CD4
model <- "~ factor(cd4cl) + factor(CHICflag) + factor(CHICflag)*factor(cd4cl)"

###- test interaction age*cd4
# model <- "~ factor(agecl) + factor(cd4cl) + factor(CHICflag) + factor(agecl)*factor(cd4cl)"

##- continuous
# model <- "~ agediag + sqrt(cd4) + agediag*sqrt(cd4)"

y <- "outdegree"
full.model = paste(y, model)
reg.sum.bs(ls = list(od), reg = lm, model = full.model) 

y = "size"
full.model = paste(y, model)
reg.sum.bs(ls = list(cluster), reg = lm, model = full.model) 

y = "binclus"
full.model = paste(y, model)
####---- logistic ----
reg.sum.bs(ls = list(cluster), reg = glm, model = full.model, family = binomial(link = "logit"))
table(cluster[[1]]$binclus)

