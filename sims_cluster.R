# rm(list=ls())
####---- include ----
detail_knitr <- TRUE
source("functions.R")

####---- lib ----
library(ape)

####---- scenario ----
scenario <- c("EqualStage0", "Baseline0")

####--- list.trees ----
list.sims <- vector("list", length(scenario))
for (s in 1:length(scenario)){
  list.sims[[s]] <- list.files('RData', full.names = TRUE, 
             path = paste('data/simulations/model0-simulate', 
                          scenario[s], sep = '') )
  names(list.sims)[s] <- scenario[s]
}
str(list.sims)

####---- list.dist ----
#### load list of dist from sims
if(TRUE){
distEqualStage0FNS <- list.files('RData', full.names=T, 
                           path = 'data/simulations/model0-simulateEqualStage0-distances')
distBaseline0FNS <- list.files('RData', full.names=T,
                              path = 'data/simulations/model0-simulateBaseline0-distances')
}

####---- ucsd clustering ----####
# ucsd_hivclust
if(FALSE){
  thresholds  <-  c(0.005, 0.015, 0.02, 0.05) # c(0.005, 0.01, 0.02, 0.05, 0.1) 
  tmax <- max(thresholds) # limit of distance considered
## function: input list of dist filenames
  ucsd <- function(ldist){
    
    for (i in 1:length(ldist)){
    ##- some processing
    load(ldist[i])
    name.sim <- substr(ldist[i], 
                       regexec(".RData",ldist[i])[[1]][1] - 5, 
                       regexec(".RData", ldist[i])[[1]][1] -1)
    folder.sim <- substr(ldist[i], regexec("-simulate", ldist[i])[[1]][1] + 9, regexec("-distances", ldist[i])[[1]][1] -1)
    dd <- as.data.frame(t(D))
    names(dd) <- c('ID1', 'ID2', 'distance')
    ##- temp el
    temp.el.fn <- paste(tempdir(), "/", name.sim, "_el.rds", sep = '')
    saveRDS(dd, file = temp.el.fn )
    ##- ucsd hivclust
    ucsd_hivclust(path.el = temp.el.fn,
                          thr = thresholds, 
                          k = tmax, 
                          out = paste("data/sim_ucsd_results/", folder.sim, '/', sep = '' ) )
    }
  } 
  
  lapply(distBaseline0FNS, ucsd)
  lapply(distEqualStage0FNS, ucsd)
}

    
####---- list.hivclust ----
###  get csv of clusters assignments in one list ###
  # getwd()
csvEqualStage0FNS <- list.files("csv", full.names=T, 
                          path = 'data/sim_ucsd_results/EqualStage0')
csvBaseline0FNS <- list.files("csv", full.names=T, 
                          path = 'data/sim_ucsd_results/Baseline0')
##- function n = 1; m = 1
list.hivclust <- function(list.csv){
  ## Structure threshold > trees
  thresholds <- c("0.005", "0.015", "0.02", "0.05") # c("0.01", "0.02", "0.05")
  ## empty list of thresholds
  cl2 <- vector("list", length(thresholds))
  ## loop
    for (m in 1:length(thresholds) ){ # index of thr
     ## vector of csv at different tree for one thr
      filenames <- list.csv[grep(thresholds[m], list.csv)]
      ## empty list of different trees
      t <- vector("list", length(filenames))
      for (n in 1:length(filenames) ){ # index of trees
        t[[n]] <- read.csv( filenames[n] )
        names(t)[n] <- substr(filenames[n], 
                              regexec("/d", filenames[n])[[1]][1] + 2,
                              regexec("_ucsd_hivclust_output", filenames[n])[[1]][1] -1)
      }
      cl2[[m]] <- t
      names(cl2)[m] <- thresholds[m]
    }
return(cl2)
  }
  
  cl_Baseline0 <- list.hivclust(list.csv = csvBaseline0FNS )
  cl_EqualStage0 <- list.hivclust(list.csv = csvEqualStage0FNS )
 # names(cl_EqualStage0[[1]]) 
  # str(cl_Baseline0[[1]][1])
  
# saveRDS(cl_Baseline0, file = "data/sim_ucsd_results/list.hivclust.sim.Baseline0.rds")
# saveRDS(cl_EqualStage0, file = "data/sim_ucsd_results/list.hivclust.sim.EqualStage0.rds")

####---- load list.hivclust ----
  cl_Baseline0 <- readRDS( file = "data/sim_ucsd_results/list.hivclust.sim.Baseline0.rds" )
  cl_EqualStage0 <- readRDS( file = "data/sim_ucsd_results/list.hivclust.sim.EqualStage0.rds" )
  
####---- heleper functions ----
  deme2age <- function(deme){ as.numeric(
    substr(regmatches( deme , regexpr( '\\.age[0-9]', deme ) ), 5, 5) 
  ) }
  deme2stage <- function(deme){as.numeric( 
    substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
  ) }
  deme2risk <- function(deme){as.numeric( 
    substr(regmatches( deme , regexpr( 'riskLevel[0-9]', deme )), 10, 10)
  ) }
  deme2care <- function(deme){as.numeric( 
    substr(regmatches( deme , regexpr( 'care[0-9]', deme )), 5, 5)
  ) }
  
  
####---- clust.stats ----
if(FALSE){
 
  ### --- start function clust.stats --- ###
  ##- calculate both numclus and sizeclus for each seqindex into a LIST
  ##- with same variable names
  ##- comprises demes states
  ##- For UCSD files, no ID if no cluster (size < 2) !
# clus=cl2[[1]]; i=1; tree = load(list.trees[["Baseline0"]][1]) ; str(clus[[1]])

  
  clust.stats <- function(clus, sim){
    
    load(sim)
    
    ## get ALL tips names and demes
    tip.names <- data.frame("id" = sub("simt_", "", daytree$tip.label),
                            "stage" = sapply( sampleDemes, deme2stage ),
                            "age" = sapply( sampleDemes, deme2age ),
                            # "care" = sapply( sampleDemes, deme2care ),
                            "risk" = sapply( sampleDemes, deme2risk ),
                            stringsAsFactors = FALSE)
    
    ## get size of clusters
    freqClust <- lapply(clus, function(x){
      as.data.frame(table(x$ClusterID),
                    stringsAsFactors = FALSE)
    })
    
    ## empty list
    ll <- list()
    
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
      # b[sample(nrow(b),10),]
      
      ##- pseudo ClusterID when size = 1
      ## starting from max(ClusterID+1)
      ## important for down-sampling
      .start <- max(b$ClusterID, na.rm = T)
      .n <- dim(b[is.na(b$ClusterID),] )[1]
      b$ClusterID[is.na(b$ClusterID)]  <- .start + 1:.n 
      
      #- colnames
      colnames(b)[which(colnames(b) =="Freq")] <- "size"
      ll[[i]] <- b
      names(ll)[i] <- names(freqClust[i])
      
    }
    return(ll)
  }
  ### --- end function clust.stats --- ###
  
  l_Baseline0 <- lapply(cl_Baseline0, 
                        function(x) clust.stats(clus = x, sim = list.sims[["Baseline0"]][1] ))
  l_EqualStage0 <- lapply(cl_EqualStage0, 
                        function(x) clust.stats(clus = x, sim = list.sims[["EqualStage0"]][1] ))
  
#   names(l_Baseline0)
#   names( l_Baseline0[[1]][[1]] )
#   length(l_Baseline0[[1]])
  
 
  ####---- saved listUKclus ----
  # saveRDS(l_Baseline0, file = "data/sim_ucsd_results/list.sim.ucsd.Baseline0.rds" )
  # saveRDS(l_EqualStage0, file = "data/sim_ucsd_results/list.sim.ucsd.EqualStage0.rds" )
  l_Baseline0 <- readRDS(file = "data/sim_ucsd_results/list.sim.ucsd.Baseline0.rds" )
  l_EqualStage0 <- readRDS(file = "data/sim_ucsd_results/list.sim.ucsd.EqualStage0.rds" )
}  


###--- stats cluster ---###
# listUKclus <- l_bs_uk ## without patients data

####---- n clusters 2 ----
## number of different clusters (counting size 1)
sapply(l_Baseline0, function(x){
  summary(sapply(x, function(x) {
    length(unique(x$ClusterID) )
  }))
})

####---- membership ----
## proportion of cluster membership
sapply(l_Baseline0, function(x){
  summary(sapply(x, function(x) sum(x$binclus) / length(x$binclus)))
})

####---- mean size 2 ----
## stats of mean size
sapply(l_Baseline0, function(x){
  summary(sapply(x, function(x) mean(x$size)))
})

####---- median size 2 ----
## stats of median size
sapply(l_Baseline0, function(x){
  summary(sapply(x, function(x) median(x$size)))
})

####---- plots ----

tab1 <- l_Baseline0[["0.05"]][[40]]
  tab0 <- l_EqualStage0[["0.05"]][[40]]
str(tab1)  
##- test cluster size 1 -> 0
# tab1[tab1$size == 1, ]$size <- NA
boxplot(tab1$size ~ tab1$age) 
boxplot(tab1$size ~ tab1$risk) 
boxplot(tab1$size ~ tab1$stage)
############# test function downsample

###- downsample in a list
downsample <- function(df, iter = 2){
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
    names(down_listclus)[i] <- i
  }
  return(down_listclus)
}
###- downsample in a df
downsample2 <- function(df, iter = 2){
  ##- sampling one id per cluster k times
  k <- iter
  ## loop
  ## empty dataframe
  down <- data.frame()
  ##  k selection of one id from df
  ##   sampled by each ClusterID
  for (i in 1:k){
    down <- rbind(down, df[df$id %in% 
                               tapply(df$id,
                                      df$ClusterID, 
                                      function(x) sample(x, 1)),]
    )
  
  }
  return(down)
}

# down_tab1 <- downsample(tab1, iter = 2)
# names(down_tab1)

down_tab2 <- downsample2(tab1, iter = 30)
boxplot( down_tab2$size ~ down_tab2$age)
boxplot( down_tab2$size ~ down_tab2$risk)
boxplot( down_tab2$size ~ down_tab2$stage)
aggregate(down_tab2$size, by = list("age" = down_tab2$age), summary)
aggregate(down_tab2$size, by = list("risk" = down_tab2$risk), summary)
aggregate(down_tab2$size, by = list("stage" = down_tab2$stage), summary)

u.test <- function(df){
  U <- wilcox.test(df[df$risk == 2, "size"], 
            df[df$risk == 1, "size"], 
            alternative = "greater") # "two.sided", "less"
  return(U$p.value)
}
mean(sapply(l_Baseline0[["0.015"]], u.test) < 0.05)
mean(sapply(l_EqualStage0[["0.015"]], u.test) < 0.05)

## after downsampling

down_baseline <- lapply(l_Baseline0[["0.015"]], function(x) downsample(x, iter = 10))
down_equal <- lapply(l_EqualStage0[["0.015"]], function(x) downsample(x, iter = 10))

names(down_baseline[[1]])
mean(sapply(down_baseline, function(x) sapply(x, u.test)) < 0.05)
mean(sapply(down_equal, function(x) sapply(x, u.test)) < 0.05)

# head(down_listclus[[1]])
  # dim(down_listclus[[1]])
 
####---- run down-sample 1 ----
### over thresholds

## ordinal variables
dd_ord <- sapply( listclus, function(x) {
  downsample(df = x, iter = 100)
})


###############
    
####---- individual bootstrap regressions ----
###--- start function
###- summarize regression on bootstrap
reg.sum.bs <- function(ls, reg, model, alpha = 0.05, ...){
  
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
  }
  )
  
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
    
    return(list("model" = model, "mean parameter" = mean.parms, "signif pvalue" = sum.signif, "mean r.squared" = mean.r2)) 
  } else {
    
    return(list("model" = model, "mean parameter" = mean.parms, "signif pvalue" = sum.signif))
  }
}
###--- end function 

model1 <- "scale(size) ~ factor(stage)"
model2 <- "scale(size) ~ factor(age)"
model3 <- "scale(size) ~ factor(risk)"
model4 <- "scale(size) ~ factor(age) + factor(stage) + factor(age)*factor(stage)"
model5 <- "scale(size) ~ factor(risk) + factor(stage) + factor(risk)*factor(stage)"
models <- as.list(paste0("model", 1:5))

## example
## age
reg.sum.bs(ls = l_Baseline0, reg = lm, model = model2) 
reg.sum.bs(ls = l_EqualStage0, reg = lm, model = model2) 
## age and stage
reg.sum.bs(ls = l_Baseline0, reg = lm, model = model4) 
reg.sum.bs(ls = l_EqualStage0, reg = lm, model = model4) 

####---- all lm ----
lapply(models, function(x) {reg.sum.bs(ls = l_Baseline0, reg = lm, model = x)
  })
lapply(models, function(x) {reg.sum.bs(ls = l_EqualStage0, reg = lm, model = x)
})

####---- logistic ----
reg.sum.bs(ls = listUKclus, reg = glm, model = logit_model_uk, family = binomial(link = "logit"))

####---- stop ----

## TODO
## revise downsampling with mean (or median) of covariates by cluster (and not cluster size)
## add thresholds: change ucsd function to do nothing if csv files exist at all thresholds
## categorize CD4 and age as for SA
## 

# ADD LSD results
# ADD SA method results

names(l_Baseline0)
test1 <- l_Baseline0[["0.015"]][[1]]
test0 <- l_EqualStage0[["0.015"]][[1]]
str(test1)

by(test1$size, test1$risk, mean)
by(test0$size, test0$risk, mean)

boxplot(test1$size ~ test1$risk)
boxplot(test0$size ~ test0$risk)

### without downsample
sizes_stage1 <-  test[test$stage == 1, ]$size 
sizes_otherstages <-  test[test$stage != 1, ]$size 
dev.off()
plot (density(sizes_stage1))
lines (density(sizes_otherstages))
wilcox.test(sizes_stage1, sizes_otherstages, alternative = 'greater')

## with downsample
head(test1)
m <- aggregate(test1[, c("stage", "age", "risk")], 
               by = list("ClusterID" = test1$ClusterID, "size" = test1$size), 
                    FUN = mean)

str(m)
tail(m)
plot(m$stage, m$size)
plot(m$age, m$size)
plot(m$risk, m$size)

# U test for the equal rates simulations: 
wtest_er <- lapply( l_Baseline0, function(obs) {
  wilcox.test( obs[[1]] / pstage[1] # NOTE pstage[1] is propto average duration of EHI
               , obs[[5]] / (1-pstage[1] )  #NOTE 1-pstage[1] is propto average duration of the rest of the infectious period
               , alternative = 'greater' #NOTE this is a one-tailed test. H1: EHI rate > !EHI rate
  )
})

# U test for the baseline simulations: 
wtest_bl <- lapply( obs_bl, function(obs) {
  wilcox.test( obs[[1]] / pstage[1] 
               , obs[[5]] / (1-pstage[1] ) 
               , alternative = 'greater'
  )
})
##########################################
##########---- laaaaa  ----###############
##########################################
##########################################
