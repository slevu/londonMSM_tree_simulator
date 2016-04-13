#########################################################
###--- use SA method from phydynR on UK bootstrap trees
###--- dated with LSD
#########################################################

# rm(list=ls())
library(ape)
library(phydynR)
source("functions.R")

 ##- list of lsd trees
 ## filename pattern from LSD changes with LSD version !
 # list.lsd.trees <- list.files(path = "data/LSD", pattern = "result.date", full.names = TRUE) # macbook
   list.lsd.trees <- list.files(path = "data/LSD", pattern = "result_newick_date", full.names = TRUE)
 head(list.lsd.trees)
 
 ##- read first LSD tree to name 
 ## sampling times and CD4s with tip.labels
 t <- read.tree( list.lsd.trees[2] )
 str(t)
 STFN <- "data/LSD/t000.dates"

 ##- parameter phydynR
 ##- Maximum height
 MH <- 20 # look up to 10 years in past for infector probs
 ##- CD4 values
 load("../phylo-uk/data/sub.RData")
 rm(s)
   ## selection of df covariates
   names(df)
   cd4s <- setNames(df$cd4, df$seqindex)[t$tip.label]
   head(cd4s)
 ##- incidence, prevalence: central scenario # todo: range of values
 ## Yin et al. 2014: 2,820 (95% CrI 1,660-4,780)
 newinf <- 2500 # c(1660, 4780)
 plwhiv <- 43150 / 2 # c(43510 / 2, 43510 / 1.5) 
 ##- sampling times
 dates <-( read.table(STFN, skip=1, 
                      colClasses=c('character', 'numeric') ) )
 head(dates)
 ##- named vector
 sampleTimes <- setNames( dates[,2], dates[,1] )[t$tip.label] 
 head(sampleTimes)
 
###--- function to apply SA over lsd tree
sa <- function(lsd_tree){
  W <- phylo.source.attribution.hiv( lsd_tree, 
          sampleTimes, # must use years
          cd4s = cd4s, # named numeric vector, cd4 at time of sampling
          ehi = NA, # named logical vector, may be NA, TRUE if patient sampl
          numberPeopleLivingWithHIV = plwhiv, # scalar
          numberNewInfectionsPerYear = newinf, # scalar
          maxHeight = MH,
          res = 1e3,
          treeErrorTol = Inf)
  return(W)
}

## loop
for (i in 1:length(list.lsd.trees)){
 tree <- read.tree(file = list.lsd.trees[i])
 W <- sa(lsd_tree = tree)
 saveRDS(W, file = paste("data/phydynR/W0_uk_mh", MH, "_",  i, ".rds", sep = '') )
}

####---- reprise ----####
###--- Do these things to outdegree and cluster sizes :
## add explanatory variables
## categorize age and cd4
## run model and tests

## list of infector prob files
list.W0 <- list.files("data/phydynR", pattern = 'mh20', full.names = TRUE)
## order
list.W0 <- list.W0[order(nchar(list.W0), list.W0)]
 
###--- analyses ---###
##- calculating outdegree

### list of outdegrees for m bootstrap
### function: input filename of W
outdegree <- function(w.fn){
  W <- readRDS(w.fn)
  out <- aggregate(
    x = list(outdegree = W$infectorProbability),
    by = list(patient = W$donor), 
    FUN = function(x) sum(x, na.rm = T) )
  return(out)
} 

list.outdegree <- lapply(list.W0, outdegree)
## nb patients
summary(sapply(list.outdegree, function(x) dim(x)[1] ))
## head
head(list.outdegree[[1]])
## summary OD
summary(list.outdegree[[1]]$outdegree)

###--- clusters ---###
## get cluster size for one threshold (say 2%) in parallel with outdegree before doing regressions
l_bs_uk <- readRDS( file = "data/listUK_ucsd_clus.rds")
## thr
thr <- "0.02"
## list of m bootstrap dataframe of cluster assignements
list.clus <- l_bs_uk[[thr]]
str(list.clus[[1]])

### add demo variables ...
##- add individual explanatory variates
##- selection of df covariates
y <- df[,c("seqindex","patientindex", 
           "dob_y", "agediag", "cd4", 
           "ydiag", "CHICflag", "ethnicityid")]

y$ethn.bin <- ifelse(y$ethnicityid == "White", "white", "not white")
y$CHICflag <- ifelse(y$CHICflag == "Yes", 1, 0)
y$ethnicityid <- NULL
y <- unfactorDataFrame(y)
str(y)

## categorize continuous variables
# str(out)
age2quantile <- function(age){
  if (is.na(age)) return (NA)
  if (age < 27) return(1)
  if (age < 33) return(2)
  if (age < 40) return(3)
  return(4)
}

cd4toStage <- function(cd4){
  if (is.na(cd4)) return(NA)
  if (cd4 > 700 ) return(1) #based on .9 quantile
  if (cd4 > 500 ) return(2)
  if (cd4 > 350 ) return(3)
  if (cd4 > 200 ) return(4)
  return(5)
}

y$agecl <- sapply( y[ , "agediag"] , age2quantile )
y$cd4cl <- sapply( y[ , "cd4"] , cd4toStage )
head(y)

### add demo variables
cluster <- lapply(list.clus, function(x){
  merge(x, y, 
        by.x = "id", by.y = "seqindex", 
        all.x = T, sort = FALSE)
  })

od <- lapply(list.outdegree, function(x){
  merge(x, y, 
        by.x = "patient", by.y = "seqindex", 
        all.x = T, sort = FALSE)
})

##- function to prune cluster according to a threshold of sampling time to control for cohort effect ?
## or use ydiag ? 
## depends on clustering algorithm ?
a <- cluster[[1]][, c("id", "ClusterID", "ydiag")]
prune.clus <- function(a, thr = 1900){
  ## subset df by year
  a <- a[a$ydiag > thr,]
  ## for each clusterID, assign new size
  for (i in unique(a[,"ClusterID"])){
    a[ a[,"ClusterID"] == i, "size" ] <- nrow(a[ a[,"ClusterID"] == i, ])
  }
  return(a)
}

b <- prune.clus(a, 1900)

summary(b$size)
summary(cluster[[1]]$size)
plot(b$size[order(b$size)],
          cluster[[1]]$size[order(cluster[[1]]$size)])

## plot on first bootstrap
## 
boxplot(sqrt(od[[1]]$outdegree) ~ od[[1]]$agecl)
boxplot(sqrt(od[[1]]$outdegree) ~ od[[1]]$cd4cl)
boxplot(sqrt(od[[1]]$outdegree) ~ od[[1]]$CHICflag)

boxplot(sqrt(cluster[[1]]$size) ~ cluster[[1]]$agecl)
boxplot(sqrt(cluster[[1]]$size) ~ cluster[[1]]$cd4cl)
boxplot(sqrt(cluster[[1]]$size) ~ cluster[[1]]$CHICflag)

############################################
###########-------- la ------############### 
############################################    

##- downsample ?

reg.sum.bs
reg.sum

summary(lm("size ~ scale(CHICflag) + scale(sqrt(cd4)) ", cluster[[1]]))

summary(
  glm("binclus ~ factor(CHICflag) ",
      family = binomial(link = "logit"), cluster[[3]])
  )
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

###--- end ---###