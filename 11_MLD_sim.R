##---- MLD regression on simulations ----
##- find estimates of W_ij
##- define incident infection
##- find MLDs: sample with infector probabilities
##- correlates of being MLD vs not (after having tested zero-inflated regressions)
##- on two simulation scenarios
##- summarize distribution of coefficients and p-values

##- for cluster:
##- MLD sampled from neighborhood

# rm(list=ls())
##---- libs ----
library(pscl) # zeroinfl
library(MASS) # glm.nb
library(lmtest)
library(Hmisc)

##---- sim files ----
source("load_sims_files.R")

##---- helper functions ----
{
  deme2age <- function(deme){ as.numeric(
    substr(regmatches( deme , regexpr( '\\.age[0-9]', deme ) ), 5, 5) 
  ) }
  deme2stage <- function(deme){as.numeric( 
    substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
  ) }
  deme2risk <- function(deme){as.numeric( 
    substr(regmatches( deme , regexpr( 'riskLevel[0-9]', deme )), 10, 10)
  ) }
  ##- remove outliers
  winsorize2 <- function (x, multiple=3)
  {
    if(length(multiple) != 1 || multiple <= 0) {
      stop("bad value for 'multiple'")
    }
    med <- median(x)
    y <- x - med
    sc <- mad(y, center=0) * multiple
    y[ y > sc ] <- sc
    y[ y < -sc ] <- -sc
    y + med
  }
  ##-
  unfactorDataFrame <- function( x ) {
    x <- data.frame( lapply(x, as.character), stringsAsFactors = FALSE)
    x <- data.frame( lapply(x, type.convert, as.is = TRUE),
                     stringsAsFactors = FALSE)
    return(x)
  }
} ## functions

##---- check id ----
## distribution of in-degree for incident
check_id <- function() {
  load(list.sims[['Baseline0']][1])
  ## get ALL tips names and demes
  tip.states <- data.frame(
    "id" = bdt$tip.label,
    "stage" = sapply( sampleDemes, deme2stage ),
    "age" = sapply( sampleDemes, deme2age ),
    "risk" = sapply( sampleDemes, deme2risk ),
    stringsAsFactors = FALSE)
  rownames(tip.states) <- NULL
  
  ##- select ID of stage 1
  index_incident <- tip.states[which(tip.states$stage == 1), 'id']
  W <- as.data.frame(W, stringsAsFactors = FALSE)
  names(W)[which(names(W) == 'infectorProbability')] <- 'ip'
  Winc <- W[W$recip  %in% index_incident, ]
  
  id <- tapply(W$ip, W$recip, sum)
  print('in-degrees for all recipients')
  print(summary(id))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.0000  0.1913  0.4458  0.4422  0.6773  1.1770 
  id_inc <- tapply(Winc$ip, Winc$recip, sum)
  print('in-degrees for stage-1 recipients')
  print(summary(id_inc)) #
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.00000 0.01153 0.16030 0.25560 0.44600 0.96200 
  
  ##- plot
  plot (density(id_inc), col = 'red', main = 'In-degrees', xlab = '')
  lines (density(id), col = 'blue')
  abline(v = mean(id_inc), col = "red", lty = 2)
  abline(v = median(id), col = "blue", lty = 3)
  legend("topright", legend = c('Recipients in stage 1', 'All recipients'), fill = c('red', 'blue'))
}  ##- compute and plot in-degree mean, median
check_id()

##---- compute mld counts ----
## save list of tip.states and counts of MLD per tip
get.mc.mld.counts <- function(x, Nsim = 100){
  ##- test on first: load(l[1])
  load(x) # RData list of simulation results (including dated tree, sampleDemes and W_ij, ... distance, cd4s, newinf, plwhiv)

  ## get ALL tips names and demes
  tip.states <- data.frame(
    "id" = bdt$tip.label,
    "stage" = sapply( sampleDemes, deme2stage ),
    "age" = sapply( sampleDemes, deme2age ),
    "risk" = sapply( sampleDemes, deme2risk ),
    stringsAsFactors = FALSE)
  rownames(tip.states) <- NULL
  
  ##- select ID of stage 1 as recipients and corresponding infector probs
  index_incident <- tip.states[which(tip.states$stage == 1), 'id']
  W <- as.data.frame(W, stringsAsFactors = FALSE)
  names(W)[which(names(W) == 'infectorProbability')] <- 'ip'
  Winc <- W[W$recip  %in% index_incident, ]
  
  source("single.mc.mld.R")
  ##- apply MLD Nsim times
  mc.mld <- replicate(Nsim, single.mc.mld.ip(x = Winc))
  
  ##- get rid of row = donor 'out'
  counts.mld <- mc.mld[-which(rownames(mc.mld)=='out'), ]
  return(list('states' = tip.states, 'counts' = counts.mld))
}

##---- run MC ----
MC.COUNT_ba = paste(path.results, 'mc.mld.baseline.rds', sep = '/') #"./data/mc.mld.baseline.rds"
MC.COUNT_er =  paste(path.results, 'mc.mld.equalrates.rds', sep = '/') #"./data/mc.mld.equalrates.rds"
mc <- function(LST, FN){
  if(!file.exists(FN)){
    s <- system.time(
      list.counts.mld <- lapply(LST, get.mc.mld.counts)
    ) # 1979s = 33mn for 100 sims * 100 MC reps
    print(s)
    saveRDS(list.counts.mld, FN)
  } else {
    list.counts.mld <- readRDS(FN)
  }
  return(list.counts.mld)
}
list.counts.mld_ba <- mc(LST = list.sims[['Baseline0']], FN = MC.COUNT_ba)
list.counts.mld_er <- mc(LST = list.sims[['EqualStage0']], FN = MC.COUNT_er)

##---- stats
stats <- function(list.counts.mld){ # stats
  PropOfNonZero <- lapply(list.counts.mld, function(x) round(apply(x$counts, 2, function(z) sum(z>0)) / dim(x$counts)[1] * 100, 2))
  s <- sapply(PropOfNonZero, summary) # 4% on non-zero
  m <- sapply(list.counts.mld, function(x) max(x$counts))
  
  # ZZ <- lapply(list.counts.mld, function(x) apply(x$counts, 2, function(z) setNames(as.vector(table(z)), names(table(z)))))
  # ZZZ <- lapply(ZZ, function(x) tapply(do.call(c,x), names(do.call(c,x)), sum))
  # 
  # mm <- do.call(rbind, lapply(list.counts.mld, function(x) x['counts']))
  # mmt <- lapply(mm[1:10], function(x) as.vector(prop.table(table(x))))
  list(prop_non_zero = summary(s['Mean',]), stats_max = summary(m))
} ## stats counts
sb <- stats(list.counts.mld = list.counts.mld_ba)
se <- stats(list.counts.mld = list.counts.mld_er)

##---- test ----
if(FALSE){
  ##-- 10 reps
  ##- test: bind 10 random MC replicates
  ## debug: 
  counts.mld <- list.counts.mld_ba[[1]][['counts']]; tip.states <- list.counts.mld_ba[[1]][['states']]
  ks <- sample(1:100, 10) # 1:100
  tab.mld.10 <- data.frame('id' = rep(rownames(counts.mld), length(ks)), 'mld' = as.vector(counts.mld[, ks]))
  table(tab.mld.10$mld)
  tab <- merge(x = tip.states, 
               y = tab.mld.10, by = 'id', all.x = T) # all.x = T)
  table(tab$mld, useNA = 'ifany')
  tab$mld[is.na(tab$mld)] <- 0L # all donors (if all.x = T)
  ##- test: bind 10 random MC replicates
  
  ##-- plots 
  par(mfrow = c(1, 3))
  boxplot(stage ~ mld, data = tab, horizontal =TRUE, ylab = 'MLD count', xlab = 'stage')
  boxplot(age ~ mld, data = tab, horizontal = T, xlab= 'age sampling')
  boxplot(risk ~ mld, data=tab, horizontal = T, xlab = 'risk level')
} # plot on 10 reps

##- question: do zero-inflated poisson or neg binomial models provide better fit than logistic ?
##- answer: no
##---- compare reg ----
{
  ##---- first replicate ----
  ##- keep one df, merge copy MLD count
  counts.mld <- list.counts.mld_ba[[1]][['counts']]; tip.states <- list.counts.mld_ba[[1]][['states']]
  k <- 1
  tab.mld <- data.frame('id' = rownames(counts.mld), 'mld' = counts.mld[, k])
  tab <- merge(x = tip.states,
               y = tab.mld, by = 'id', all.x = T) # all.y = T)
  tab$mld[is.na(tab$mld)] <- 0L # all donors (if all.x = T)
  tab$binmld <- ifelse(tab$mld == 0, 0, 1)
  tab$ehi <- ifelse(tab$stage == 1, 1, 0)
  tab$young <- ifelse(tab$age == 1, 1, 0)
  #with( tab, table(mld, age, stage) ) # 'in the general realm of very rare events'
  #with( tab, table(mld, ehi, young) ) # 'in the general realm of very rare events'
  
  ##---- models ----
  # X4 <- c('factor(age)', 'factor(stage)', 'risk')
  X4 <- c('young', 'ehi', 'risk')
  Y = 'mld'
  Z = 1
  mod1 <- as.formula(paste(Y, '~', paste(X4, collapse = ' + ')))
  mod1bi <- as.formula(paste('binmld', '~', paste(X4, collapse = ' + ')))
  
  ##---- fit ----
  logit1 <- glm(mod1bi, family="binomial", data = tab)
  poisson1 <- glm(mod1 , family="poisson", data = tab) # poisson GLM
  nb1 <- glm.nb(mod1, data = tab) # negative binomial
  zi_poisson1 <- zeroinfl(mod1, data = tab, dist= 'poisson') # zero-inflated poisson
  tryCatch(zi_nb1 <- zeroinfl(mod1, data = tab, dist= 'negbin'), error = function(e) e) # zero-inflated neg bin
  
  ## coefs
  fm <- list("LOGIT" = logit1, "PO" = poisson1, "NB" = nb1, "ZIPO" = zi_poisson1, "ZINB" = zi_nb1)
  k <- which.max( sapply(fm, function(x) length(coef(x))) )
  full.names <- names(coef(fm[[k]])) 
  coefs <- sapply(fm, function(x) coef(x)[1:length(full.names)])
  rownames(coefs) <- full.names
  ## for latex
  tt <- data.frame(coefs, check.names = FALSE)
  
  
  # pseudo-code for extracting p-values in both configuraiton of coef
  get.glm.pvalues <- function(x) {
    if (any(class(x) == "glm")){
      p <- summary(x)$coef[,4]
      names(p) <- paste0('count_', names(p))
    } else {
      l <- lapply(summary(x)$coef, function(x) x[,4])
      p <- lapply(1:2, function(d) setNames(l[[d]], paste0(names(l)[d], '_', names(l[[d]]))) )
      p <- unlist(p)
    }
    return(p)
  }
  
  get.all.glms.pvalues <- function(list.of.fit = fm){
    k <- which.max( sapply(list.of.fit, function(x) length(coef(x))) )
    full.names <- names(coef(fm[[k]])) 
    m <- sapply(list.of.fit, function(x) get.glm.pvalues(x)[1:length(full.names)])
    rownames(m) <- full.names
    return(m)
  }
  ## p-values < 0.05
  p.sign <- get.all.glms.pvalues(fm) < 0.05
  
  ## format for latex table
  cell.fmt <- matrix(rep("", nrow(tt) * ncol(tt)), nrow = nrow(tt))
  cell.fmt[p.sign] <- "bfseries"
  #cell.fmt[p.sign] <- "color{red}"
  
  ## add likelihood
  ## log-L
  gofs <- rbind('log-likelihood' = sapply(fm, function(x) logLik(x)),
                'df' = sapply(fm, function(x) attr(logLik(x), "df")),
                AIC = sapply(fm, function(x) AIC(x)))
  tt <- rbind(tt, gofs )
  cell.fmt <- rbind(cell.fmt, matrix(rep("", nrow(gofs) * ncol(gofs)), nrow = nrow(gofs)))
  # dim(tt); dim(cell.fmt)
  nr <- dim(coefs)[1]/2 # size of rgoup
  tab.latex <- latex(round(tt,2), title = '', n.rgroup = c(nr, nr, 3), cellTexCmds = cell.fmt, booktabs = TRUE, numeric.dollar = FALSE, rowname = latexTranslate(rownames(tt)))
  
  #tab.latex
  html(tab.latex, rmarkdown = TRUE)
} ## tests zero-inflated regressions

##---- logit ----
## debug: s = 1; k = 1; counts = list.counts.mld[[s]][[2]]; df = list.counts.mld[[s]][[1]]; X =  c('factor(age)', 'factor(stage)', 'risk')
## debug:  k = 1; counts = counts.mld; df = tip.states; X =  c('factor(age)', 'factor(stage)', 'risk')
get.mld.logit <- function(k, counts = counts.mld, df = tip.states, X =  c('young', 'ehi', 'risk')){ # index column, states with ID variable
  tab.mld <- data.frame('id' = rownames(counts), 'mld' = counts[, k], stringsAsFactors = FALSE)
  tab <- merge(x = df,
               y = tab.mld, by = 'id', all.x = T) # all.y = T)
  tab$mld[is.na(tab$mld)] <- 0L # all donors (if all.x = T)
  
  ##- 0. logistic regression (0, non 0)
  mod1bi <- as.formula(paste('binmld', '~', paste(X, collapse = ' + ')))
  tab$binmld <- ifelse(tab$mld == 0, 0, 1)
  tab$ehi <- ifelse(tab$stage == 1, 1, 0)
  tab$young <- ifelse(tab$age == 1, 1, 0)
  logit1 <- glm(mod1bi, family="binomial", data = tab)
  test <- list('coef' = summary(logit1)$coef, 'Loglikelihood' = logLik(logit1), 'AIC' = AIC(logit1))
  # or <- exp(coef(logit1)) # CI too slow:   or <- exp(cbind(OR = coef(logit1), confint(logit1)))
  return(test)
}
# fits.logit <- lapply(1:dim(counts.mld)[2], function(z) get.mld.logit(k = z))

##---- run logistic ----
FITLOGIT_ba <- paste(path.results, 'fit.logit.mld.baseline2.rds', sep = '/') # './data/fit.logit.mld.baseline.rds'
FITLOGIT_er <- paste(path.results, 'fit.logit.mld.equalrates2.rds', sep = '/') # './data/fit.logit.mld.equalrates.rds'
fit.glm <- function(LST, FN){
  if(!file.exists(FN)){
    st <- system.time(
      fits.logit <- lapply(LST, function(x){
        lapply(1:dim(x$counts)[2], function(z) get.mld.logit(k = z, counts = x$counts, df = x$states)) 
      })
    ) # long
    print(st)
    saveRDS(fits.logit, FN)
  } else {
    fits.logit <- readRDS(FN)
  }
  return(fits.logit)
}
fits.logit_ba <- fit.glm(LST = list.counts.mld_ba, FN = FITLOGIT_ba)
fits.logit_er <- fit.glm(LST = list.counts.mld_er, FN = FITLOGIT_er)

##- get distribution of coefs (direction): 100 tree sims * 100 MC sims
## fit = fits.logit_ba[[1]]; fit = fits.logit_ba[31:40]
get.coef <- function(fit){
  # source('functions.R') # for winsorize2
  # m <- sapply(fit, function(x) x$coef[2:4,1]) # one sim * 100 MC sims
  m <-  do.call(cbind, lapply(fit, function(l) sapply(l, function(x) x$coef[2:4,1]) ) )
  m <- t(apply(m, 1, winsorize2)) ##- get rid of outliers
  
  ##- plot
  mxy <- sapply(apply(m, 1, density), 
                function(d) c(minx = min(d$x), maxx = max(d$x), 
                              miny = min(d$y), maxy = max(d$y))) ## find x and y ranges
  plot (density(m[1,]), xlab = '', main = 'Regression coefficients', col = 'red', 
        xlim = c(min(mxy['minx',]), max(mxy['maxx',])),
        ylim = c(min(mxy['miny',]), max(mxy['maxy',])))
  lines (density(m[2,]), col = 'blue')
  lines (density(m[3,]), col = 'green')
  legend("topleft", legend = rownames(m), fill = c('red', 'blue', 'green'))
}

##- look at distribution of coef estimate with SA method
get.coef(fits.logit_ba)
get.coef(fits.logit_er)

##- count any p-value < 0.05 by X
## debug: X = X; fit = fits.logit[[1]][[1]]
get.p.values <- function(fit, X =  c('young', 'ehi', 'risk')){ #c('factor(age)', 'factor(stage)', 'risk')){
  #m <- summary(fit)
  p.values <- fit$coef[,4] 
  ##- get rid of parentheses
  unparenthesis <- function(string){
    gsub("\\(|\\)", "", string)
  }
  ##- return list of X variables with p-value
  P <- sapply(X,  function(v) any(p.values[grep(unparenthesis(v),
                                                unparenthesis(names(p.values)) )] < 0.05)
  )
  return(P)
}

##---- stats p-values ----
## by sim
p.sum_ba <- lapply(fits.logit_ba, function(x) sapply(x, function(z) get.p.values(fit = z)))
p.sum_er <- lapply(fits.logit_er, function(x) sapply(x, function(z) get.p.values(fit = z)))
## by MC replicate
p.sum.per.sim_ba <- sapply(p.sum_ba, rowSums)
p.sum.per.sim_er <- sapply(p.sum_er, rowSums)

apply(p.sum.per.sim_ba, 1, summary)
#         factor(age) factor(stage)   risk
# Min.           1.00           100  10.00
# 1st Qu.        7.00           100  79.50
# Median        16.00           100  94.00
# Mean          23.86           100  85.57
# 3rd Qu.       36.00           100  98.00
# Max.          81.00           100 100.00
apply(p.sum.per.sim_er, 1, summary)
#         factor(age) factor(stage)   risk
# Min.           0.00         99.00  94.00
# 1st Qu.        5.00        100.00 100.00
# Median        11.00        100.00 100.00
# Mean          16.64         99.91  99.84
# 3rd Qu.       20.25        100.00 100.00
# Max.          81.00        100.00 100.00

##---***---***---***---***---***---
##---- cluster ----
# cw_Baseline0 <- readRDS(file = paste(path.results, 'list.sim.clus-outdeg.Baseline0.rds', sep = '/') )
# ##- first elements of second list level
# struc <- lapply(cw_Baseline0, `[`, 1) # or lapply(cw_Baseline0, function(x) x[1])
# str(struc[1:2], give.attr=T, give.length = F, give.head=T, vec.len = 1, indent.str="|", comp.str="----")

#- load D
#- format edgelist
#- for each incident patient,
#- select pair with minimum distance 
#- and consider other node as donor (see ugly loop)
#- OR
#- for each incident, draw one donor in same cluster and add MLD count


get.mld.clust <- function(sim, Nsim = 100) {
  load(sim) # load( list.sims[['Baseline0']][2] )
  ## get ALL tips names and demes
  tip.states <- data.frame(
    "id" = bdt$tip.label,
    "stage" = sapply( sampleDemes, deme2stage ),
    "age" = sapply( sampleDemes, deme2age ),
    "risk" = sapply( sampleDemes, deme2risk ),
    stringsAsFactors = FALSE)
  rownames(tip.states) <- NULL
  
  ##- select ID of stage 1
  index_incident <- tip.states[which(tip.states$stage == 1), 'id']
  D <- distance[['D']]
  D[,1:2] <- apply(D[,1:2], 2, as.character)
  D_inc <- D[D[,1]  %in% index_incident | D[,1]  %in% index_incident, ]
  
  source("single.mc.mld.R")
      # draw 1 MLD in neighborhood
    # counts.mld <- replicate(Nsim, single.mc.mld.clust(edgelist = D_inc, ids = index_incident, thr = 0.015))
    counts.mld <- sapply(1:Nsim, function(x) {
      print(x)
      single.mc.mld.clust(edgelist = D_inc, ids = index_incident, thr = 0.015)
    })
    return(list('states' = tip.states, 'counts' = counts.mld))
}

##---- run MC ----
MC.COUNT.CLUS_ba = paste(path.results, 'mc.mld.clus.baseline.rds', sep = '/') 
MC.COUNT.CLUS_er =  paste(path.results, 'mc.mld.clus.equalrates.rds', sep = '/') 
mc <- function(LST, FN){
  if(!file.exists(FN)){
    s <- system.time(
      list.counts.mld <- lapply(LST, get.mld.clust)
    ) # 1979s = 33mn for 100 sims * 100 MC reps
    print(s)
    saveRDS(list.counts.mld, FN)
  } else {
    list.counts.mld <- readRDS(FN)
  }
  return(list.counts.mld)
}
list.counts.mld_ba <- mc(LST = list.sims[['Baseline0']], FN = MC.COUNT.CLUS_ba)
list.counts.mld_er <- mc(LST = list.sims[['EqualStage0']], FN = MC.COUNT.CLUS_er)

#### laaaaaaaa
LST = list.sims[['Baseline0']]
list.counts.mld <- lapply(LST[1:2], get.mld.clust)
table(list.counts.mld$`10222`$counts[,1])
tail(list.counts.mld$`10222`$counts[,1])

##- test logit 
LST = list.counts.mld
fits.logit <- lapply(LST, function(x){
  lapply(1:dim(x$counts)[2], function(z) get.mld.logit(k = z, counts = x$counts, df = x$states)) 
})


##- apply: 
FITLOGITCLUS_ba <- paste(path.results, 'fit.logit.mld.clust.baseline.rds', sep = '/')
FITLOGITCLUS_ba2 <- paste(path.results, 'fit.logit.mld.clust.baseline2.rds', sep = '/')
FITLOGITCLUS_er <- paste(path.results, 'fit.logit.mld.clust.equal.rds', sep = '/')
FITLOGITCLUS_er2 <- paste(path.results, 'fit.logit.mld.clust.equal2.rds', sep = '/')

apply.mld.logit.cluster <- function(lst, fn, m = 1){
  if(!file.exists(fn)){
    ltest <- lapply(lst, function(x) fit.logit.mld.clust(x, method = m))
    saveRDS(ltest, fn)
  } else {
    ltest <- readRDS(fn)
  }
  pp <- sapply(ltest, get.p.values )
  return(apply(pp, 1, sum))
}

##- tests results
apply.mld.logit.cluster(lst = list.sims[['Baseline0']],
                        fn = FITLOGITCLUS_ba,
                        m = 1)
# young   ehi  risk 
# 33   100    54 

apply.mld.logit.cluster(lst = list.sims[['Baseline0']],
                        fn = FITLOGITCLUS_ba2,
                        m = 2)
# young   ehi  risk 
# 10    77     8
### WATCHOUT direction of EHI relation !!

apply.mld.logit.cluster(lst = list.sims[['EqualStage0']],
                        fn = FITLOGITCLUS_er,
                        m = 1)
# young   ehi  risk 
# 35    98    97 

apply.mld.logit.cluster(lst = list.sims[['EqualStage0']],
                        fn = FITLOGITCLUS_er2,
                        m = 2)
# young   ehi  risk 
# 20    45    18 


###------- lllllaaaaaa
###


##---- check outdegree by stage ----
x = l_ba[[1]]
x = l_er[[1]]
test <- function(x){
  load(x)
  
  W <- as.data.frame(W, stringsAsFactors = FALSE)
  names(W)[which(names(W) == 'infectorProbability')] <- 'ip'
  od <- tapply(W$ip, W$donor, sum)
  oddf <- data.frame('id' = names(od), od)
  
  head(oddf)
  tip.states <- data.frame(
    "id" = bdt$tip.label,
    "stage" = sapply( sampleDemes, deme2stage ),
    "age" = sapply( sampleDemes, deme2age ),
    "risk" = sapply( sampleDemes, deme2risk ),
    stringsAsFactors = FALSE)
  rownames(tip.states) <- NULL
  
  df <- merge(tip.states, oddf, by = 'id', all.y = TRUE)
  boxplot(od ~ stage, df)
  #tapply(df$od, factor(df$stage), mean)
}

####--- surplus -----



if(FALSE){
  single.mc.mld
  system.time( mld <- single.mc.mld(x = Winc) ) # 0.09s
  plot(mld[-which(names(mld)=='out')]) # zoom without 'out'
  head(mld)
} ##- compute one MC draw and plot counts
