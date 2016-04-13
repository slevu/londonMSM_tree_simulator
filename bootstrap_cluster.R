# rm(list=ls())
####---- include ----
detail_knitr <- TRUE
source("functions.R")

####---- lib ----
library(ape)

####---- list.tree.path ----
#### load list of path to Examl trees
if(TRUE){
path.trees <- "../HPC/phylo-msm-uk/in_HPC_20160404/"
.a <- list.files(path.trees)
list.tree.path <- .a[grep("result", .a)] # list of tree results
## get rid of tree 000
list.tree.path <- list.tree.path[-1]
head(list.tree.path, 3)
}



###--- Clustering ---###

####---- list.of.trees ----#### 
## (minus outgroup tips)

## filenames number
num <- sprintf("%03d", seq(1, 100))
## outgroup
og <- c( paste("Ref", 1:6, sep = ''), "HXB2" )

  first <- 1
  last <- 1 #length(list.tree.path) # 100
  ## empty list 
  list.of.trees <- vector("list", last - first + 1)
  ## index j of list elements 
  j <- 1
  for (i in first:last){
    t <- read.tree(file = paste(path.trees, 
                                list.tree.path[i], sep = ''))
    list.of.trees[[j]] <- drop.tip(t, og)
    names(list.of.trees)[j] <- paste("bs_tree_", 
                                     num[i], sep="")
    j <- j + 1
  }
  
####---- make edge lists ----####
if(FALSE){ 
  ###--- Save edge lists ---###
  test <-  function(i){
    TreeToEdgeList(t = list.of.trees[[i]],
                   name.output = names(list.of.trees)[[i]],
                   rate = 1,
                   output = "data/bootstrap/", 
                   stats = FALSE,
                   plot = FALSE)
  }
  
  system.time(
    lapply(seq_along(list.of.trees), test) 
  )
  # 9 h
}


####---- ucsd clustering ----####
# ucsd_hivclust

.b <- list.files("data/bootstrap/")
lel <- .b[grep("_el", .b)] # list of edge list
# lel <- lel[15:16] # test

# i <- 4  # tree
thresholds  <-  c(0.005, 0.01, 0.02, 0.05, 0.1) 
k <- max(thresholds) # limit of distance considered

if(FALSE){
## function  
  ucsd <- function(i){
    bs_clus <- ucsd_hivclust(path.el = paste( "data/bootstrap/",
                                                  lel[i], 
                                                  sep = '' ),
                             thr = thresholds, 
                             k = k, 
                             out = "data/ucsd_results/" )
    return(bs_clus)
  }
  
## empty list
  bs_clus <- vector("list", length(lel) )
## apply  
  system.time(
    bs_clus <- lapply(seq_along(lel), ucsd) 
  ) # 25 min
}
# i=1

####---- list.hivclust ----
###  get csv of clusters assignments in one list ###
  # getwd()
  path.csv <- "data/ucsd_results/"
  list.csv <- list.files(path.csv)
  # m <- 1
  
  ## Structure threshold > trees
  thresholds <- c("0.01", "0.02", "0.05") ## for now !
  ## empty list of 3 thresholds
  cl2 <- vector("list", length(thresholds))
  ## loop
  for (m in 1:length(thresholds) ){ # index of thr
    ## vector of csv at different tree for one thr
    filenames <- list.csv[grep(thresholds[m], list.csv)]
    ## empty list of different trees
    t <- vector("list", length(filenames))
    for (n in 1:length(filenames) ){ # index of trees
      t[[n]] <- read.csv( paste(path.csv, filenames[n], sep = '') )
      names(t)[n] <- num[n]
    }
    cl2[[m]] <- t
    names(cl2)[m] <- thresholds[m]
  }
  
  names(cl2[[1]]) 
  
# saveRDS(cl2, file = "data/ucsd_results/list.hivclust.rds")

####---- load list.hivclust ----
cl2 <- readRDS( file = "data/ucsd_results/list.hivclust.rds" )

####---- stats1 ----
  ###--- stats clusters without pseudo cluster size 1 ---###
  
  ## Size (=Freq) of each cluster
  aa <- lapply(cl2, function(x){
    lapply(x,
         function(x) as.data.frame(table(x$ClusterID),
                                   stringAsFactors = FALSE) )
  })
  names(aa)

####---- n clusters ----    
  ## number of different cluster of size > 1
  sapply(aa, function(x){
    summary(sapply(x, function(x) dim(x)[1]))
  })
  
####---- mean size ----  
  ## stats of mean size
  sapply(aa, function(x){
    summary(sapply(x, function(x) mean(x$Freq)))
  })

####---- max size ----   
  ## stats of max size
  sapply(aa, function(x){
    summary(sapply(x, function(x) max(x$Freq)))
  })
  
####---- plots size ----####
  ## mean cluster size > 1
  par(mfrow = c(2,3))
  for (i in 1:length(aa) ) {
    hist(sapply(aa[[i]], function(x) mean(x$Freq)), 
         main = names(aa)[i], 
         xlab = "mean cluster size" )
  }
  ## max cluster size > 1
  for (i in 1:length(aa) ) {
    hist(sapply(aa[[i]], function(x) max(x$Freq)), 
         main = names(aa)[i], 
         xlab = "max cluster size" )
  }
  dev.off()

  
  
  ####---- clust.stats ----
if(FALSE){
 
  ### --- start function clust.stats --- ###
  ##- calculate both numclus and sizeclus for each seqindex into a LIST
  ##- with same variable names
  ##- For UCSD files, no ID if no cluster (size < 2) !
# clus=cl2[[1]]; i=1; tree = list.of.trees[[1]] ; str(clus[[1]])

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
  
  l_bs_uk <- lapply(cl2, function(x) clust.stats(clus = x, tree = list.of.trees[[1]] ))
  names(cl2)
  names( l_bs_uk[[1]][[1]] )
  length(l_bs_uk[[1]])
  ####---- saved listUKclus ----
  # saveRDS(l_bs_uk, file = "data/listUK_ucsd_clus.rds")
  l_bs_uk <- readRDS( file = "data/listUK_ucsd_clus.rds")
  
  ## convert to data.frame ???
  # http://stackoverflow.com/questions/4512465/what-is-the-most-efficient-way-to-cast-a-list-as-a-data-frame
  # http://stackoverflow.com/questions/26193028/r-converting-nested-list-to-dataframe-and-get-names-of-list-levels-als-factors
  
  ##- add demo outcome: stage of infection, treatment status, age group, CHIC or not
  load("../phylo-uk/data/sub.RData")
  rm(s)
  ##- selection of df covariates
  names(df)
  # table(df$worldregion_birth, useNA = "ifany")
  
  y <- df[,c("seqindex","patientindex",  
             "agediag", "cd4", "ydiag", "CHICflag", "onartflag",
             "status", "ethnicityid")]
  
  ### recode ethnicity as character
  # table(y$ethnicityid, useNA = "ifany")
#   y$ethn <- ifelse( grepl("Black", y$ethnicityid),
#                     "Black",
#                         ifelse(grepl("Other", y$ethnicityid)|
#                                  grepl("Indian", y$ethnicityid),
#                                "Other",
#                           as.character(y$ethnicityid)))
  y$ethn.bin <- ifelse(y$ethnicityid == "White", "white", "not white")
  
  # table(y$ethnicityid, y$ethn)
  # head(y)
   y <- y[, c(1:6, 10)]
   str(y)
   # table(y$CHICflag, useNA="ifany")
   
   listUKclus <- lapply(l_bs_uk, function(x){
    lapply(x, function(x){merge(x, y, 
                    by.x = "id", by.y = "seqindex", 
                    all.x = T, sort = FALSE)
      })
  })
   
}
  
####---- saved listUKclus ----
# saveRDS(listUKclus, file = "data/listUKclus.rds")
listUKclus <- readRDS( file = "data/listUKclus.rds")
####---- stop ----




###--- stats cluster ---###

names(listUKclus)
x <- listUKclus[[1]][[2]]
summary(x$size)
str(x)

####---- n clusters 2 ----
## number of different clusters (counting size 1)
sapply(listUKclus, function(x){
  summary(sapply(x, function(x) {
    length(unique(x$ClusterID) )
  }))
})

####---- membership ----
## proportion of cluster membership
sapply(listUKclus, function(x){
  summary(sapply(x, function(x) sum(x$binclus) / length(x$binclus)))
})

####---- mean size 2 ----
## stats of mean size
sapply(listUKclus, function(x){
  summary(sapply(x, function(x) mean(x$size)))
})

####---- median size 2 ----
## stats of median size
sapply(listUKclus, function(x){
  summary(sapply(x, function(x) median(x$size)))
})

####---- dfUKclus ----####
###--- merge bootstrap results into
### a list of threshold times big dataframe
if(FALSE){
   dfUKclus <- lapply(listUKclus, function(x){
   do.call("rbind", x)
 })
 names(dfUKclus)
 str(dfUKclus)

 ####---- regression altogether ----
 ## test on regrouped bootstrap data

 lapply(dfUKclus, function(x){
   summary(lm(model0, data = x))
 })

}

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


# lm_model_uk = "scale(size) ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn) + factor(CHICflag)"
 logit_model_uk = "binclus ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn.bin) + factor(CHICflag)"
# ## example
# coef(summary(lm(lm_model_uk, data = listUKclus[[1]][[1]])))

model0 <- "scale(size) ~ scale(agediag)"
model1 <- "scale(size) ~ scale(sqrt(cd4))"
model2 <- "scale(size) ~ factor(ethn.bin)"
model3 <- "scale(size) ~ factor(CHICflag)"
model4 <- "scale(size) ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn.bin) + factor(CHICflag)"
models <- as.list(paste0("model", 0:4))

## example
# reg.sum.bs(ls = listUKclus, reg = lm, model = lm_model_uk)


####---- all lm ----
lapply(models, function(x) {reg.sum.bs(ls = listUKclus, reg = lm, model = x)
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

##########################################
##########---- laaaaa  ----###############
##########################################