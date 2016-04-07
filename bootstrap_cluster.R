# rm(list=ls())
library(ape)
source("functions.R")

## tree path
if(TRUE){
path.trees <- "../HPC/phylo-msm-uk/in_HPC_20160404/"
a <- list.files(path.trees)
l <- a[grep("result", a)] # list of tree results
## get rid of tree 000
l <- l[-1]
}

## filenames number
num <- sprintf("%03d", seq(1, 100))

## outgroup
og <- c( paste("Ref", 1:6, sep = ''), "HXB2" )

###--- Clustering ---###

##- list of trees (minus outgroup tips)

  first <- 1
  last <- 1 # 100
  ## empty list 
  list.of.trees <- vector("list", last - first + 1)
  ## index j of list elements 
  j <- 1
  for (i in first:last){
    t <- read.tree(file = paste(path.trees, l[i], sep = ''))
    list.of.trees[[j]] <- drop.tip(t, og)
    names(list.of.trees)[j] <- paste("bs_tree_", num[i], sep="")
    j <- j + 1
  }
  names(list.of.trees)
  
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


###--- ucsd clustering ---###
# ucsd_hivclust

b <- list.files("data/bootstrap/")
lel <- b[grep("_el", b)] # list of edge list
# lel <- lel[15:16] # test

# i <- 4  # tree
thresholds  <-  c(0.01, 0.02, 0.05) 
k <- max(thresholds) # limit of distance considered

if(FALSE){
## function  
  ucsd <- function(i){
    bs_clus <- ucsd_hivclust(path.el = paste( "data/bootstrap/",
                                                  lel[i], 
                                                  sep = '' ),
                             thr = thresholds, 
                             k = k, 
                             out = "ucsd_results/")
    return(bs_clus)
  }
  
## empty list
  bs_clus <- vector("list", length(lel) )
## apply  
  system.time(
    bs_clus <- lapply(seq_along(lel), ucsd) 
  ) # 25 min
}

  ###--- bin table in one list ---###
  getwd()
  path.csv <- "data/ucsd_results/"
  list.csv <- list.files(path.csv)
  # m <- 1
  
  if(FALSE){
  ## First structure trees > threshold
  ## empty list of 100 trees
  cl <- vector("list", 100)
  for (m in 1:length(num) ){ # index of trees
    ## vector of csv at different threshold for one tree
    filenames <- list.csv[grep(num[m], list.csv)]
    ## empty list at different threshold
    t <- list()
    for (n in 1:length(thresholds) ){ # index of thr
      t[[n]] <- read.csv( paste(path.csv, filenames[n], sep = '') )
      names(t)[n] <- thresholds[n]
    }
    cl[[m]] <- t
    names(cl)[m] <- num[m]
  }
  
 #names(cl[[1]]) 
  #str(cl[[1]])
  }
  
  ## Second structure threshold > trees
  ## empty list of 3 thresholds
  cl2 <- vector("list", length(thresholds))
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
  str(cl2[[1]])
  
#   ## Size (=Freq) of each cluster for thr 1
#   aa <- lapply(cl2[[1]],
#          function(x) as.data.frame(table(x$ClusterID),
#                                    stringAsFactors = FALSE) )
#   str(aa)
  
#   ## number of different cluster of size > 1
#   summary(sapply(aa, function(x) dim(x)[1]))
#   
#   ## mean cluster size
#   hist(sapply(aa, function(x) mean(x$Freq)))
  
  ## 
  
  #### llllllaaaaa  ----------------------------------------------------------
  
  ####---- clust.stats ----
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
  
  ## convert to data.frame ???
  # http://stackoverflow.com/questions/4512465/what-is-the-most-efficient-way-to-cast-a-list-as-a-data-frame
  # http://stackoverflow.com/questions/26193028/r-converting-nested-list-to-dataframe-and-get-names-of-list-levels-als-factors
  
  
  ##- add demo outcome: stage of infection, treatment status, age group, CHIC or not
  load("../phylo-uk/data/sub.RData")
  rm(s)
  ##- selection of df covariates
  names(df)
  
  
  y <- df[,c("seqindex","patientindex",  
             "agediag", "cd4", "ydiag", "CHICflag", "onartflag",
             "status", "ethnicityid")]
  
  ### recode ethnicity as character
  # table(y$ethnicityid, useNA = "ifany")
  y$ethn <- ifelse( grepl("Black", y$ethnicityid),
                    "Black",
                        ifelse(grepl("Other", y$ethnicityid)|
                                 grepl("Indian", y$ethnicityid),
                               "Other",
                          as.character(y$ethnicityid)))
  
  # table(y$ethnicityid, y$ethn)
   y <- y[, -9]
   # table(y$CHICflag, useNA="ifany")
   
   listUKclus <- lapply(l_bs_uk, function(x){
    lapply(x, function(x){merge(x, y, 
                    by.x = "id", by.y = "seqindex", 
                    all.x = T, sort = FALSE)
      })
  })
  # head(listUKclus[[1]])
  # str(listUKclus)
saveRDS(listUKclus, file = "data/listUKclus.rds")
listUKclus <- readRDS( file = "data/listUKclus.rds")

###--- regressions


###--- start function ---###
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
###--- end function ---###

str(listUKclus[[1]][[1]])
lm_model_uk = "scale(size) ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn) + factor(CHICflag)"
logit_model_uk = "binclus ~ scale(agediag) + scale(sqrt(cd4)) + factor(ethn) + factor(CHICflag)"

## example
coef(summary(lm(lm_model_uk, data = listUKclus[[1]][[1]])))

reg.sum.bs(ls = listUKclus, reg = lm, model = lm_model_uk)
reg.sum.bs(ls = listUKclus, reg = glm, model = logit_model_uk, family = binomial(link = "logit"))

## TODO
## revise downsampling with mean (or median) of covariates by cluster (and not cluster size)
