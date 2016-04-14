#########################################################
###--- use SA method from phydynR on UK bootstrap trees
###--- dated with LSD
#########################################################

# rm(list=ls())

####---- include ----
detail_knitr <- TRUE
source("functions.R")

####---- lib ----
library(ape)
library(phydynR)

####---- list of trees ----
 ##- list of lsd trees
 ## filename pattern from LSD changes with LSD version !
 if( any(grep("MacBook", Sys.info())) ){
 list.lsd.trees <- list.files(path = "data/LSD", pattern = "result.date", full.names = TRUE)
 } else {
 list.lsd.trees <- list.files(path = "data/LSD", pattern = "result_newick_date", full.names = TRUE)
 }
   
# head(list.lsd.trees)

####---- tree and sampleTImes ----
 ##- read first LSD tree to name 
 ##- sampling times and CD4s with tip.labels
 t <- read.tree( list.lsd.trees[2] )
# str(t)
 STFN <- "data/LSD/t000.dates"

 ##- CD4 values
 load("../phylo-uk/data/sub.RData")
 rm(s)
 ## selection of df covariates
# names(df)
 cd4s <- setNames(df$cd4, df$seqindex)[t$tip.label]
 head(cd4s)
 
 ##- sampling times
 dates <-( read.table(STFN, skip=1, 
                      colClasses=c('character', 'numeric') ) )
# head(dates)
 ##- named vector
 sampleTimes <- setNames( dates[,2], dates[,1] )[t$tip.label] 
 head(sampleTimes)
 
####---- params phydynR ----
 ##- Maximum height
 MH <- 20
 ##- incidence, prevalence: central scenario # todo: range of values
 ## Yin et al. 2014: 2,820 (95% CrI 1,660-4,780)
 newinf <- 2500 # c(1660, 4780)
 plwhiv <- 43150 / 2 # c(43510 / 2, 43510 / 1.5) 

####---- SA function ----
sa <- function(lsd_tree){
  W <- phylo.source.attribution.hiv( lsd_tree, 
          sampleTimes, # years
          cd4s = cd4s, 
          ehi = NA, 
          numberPeopleLivingWithHIV = plwhiv, 
          numberNewInfectionsPerYear = newinf, 
          maxHeight = MH,
          res = 1e3,
          treeErrorTol = Inf)
  return(W)
}

####---- loop SA ----
for (i in 1:length(list.lsd.trees)){
  w.fn <- paste("data/phydynR/W0_uk_mh", MH, "_",  i, ".rds", sep = '')
  if(!file.exists(w.fn)){
    tree <- read.tree(file = list.lsd.trees[i])
    W <- sa(lsd_tree = tree)
    saveRDS(W, file = w.fn )
  }
}

####---- end ----