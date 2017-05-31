##- MLD regression on simulations
##- find estimates of W_ij
##- define incident infection
##- find MLDs
##- correlates of being MLS vs not

# rm(list=ls())
##---- libs ----
library(pscl) # zeroinfl
library(MASS) # glm.nb
library(lmtest)
library(Hmisc)

##- sim files
PATHW <- "./data/simulations2/model0-simulateBaseline0"
l <- list.files(PATHW, full.names = TRUE)

##---- helper functions ----
deme2age <- function(deme){ as.numeric(
  substr(regmatches( deme , regexpr( '\\.age[0-9]', deme ) ), 5, 5) 
) }
deme2stage <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
) }
deme2risk <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'riskLevel[0-9]', deme )), 10, 10)
) }

##---- compute mld counts ----
## save list of tip.states and counts of MLD per tip
get.mc.mld.counts <- function(x, Nsim = 100){
  ##- test on first
  #load(l[1])
  load(x) # RData list of simulation results (including daytree, sampleDemes and W_ij)

  ## get ALL tips names and demes
  tip.states <- data.frame(
    "id" = daytree$tip.label,
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
  ## any NAs ?: apply(Winc, 2, function(x) sum(is.na(x)))
  
  ##- distribution of in-degree for incident
  if (FALSE){
    id <- tapply(Winc$ip, Winc$recip, sum) # id <- tapply(W$ip, W$recip, sum)
    summary(id) # mean only ~ 10%
    plot(density(id), main = 'in-degree', xlab = '', ylab = '')
    abline(v = mean(id), col = "blue", lwd = 1)
    abline(v = median(id), col = "red", lwd = 1)
  } ##- compute and plot in-degree mean, median
  
  source("../phylo-uk/code/single.mc.mld.R")
  if(FALSE){
    single.mc.mld
    system.time( mld <- single.mc.mld(x = Winc) ) # 0.09s
    plot(mld[-which(names(mld)=='out')]) # zoom without 'out'
    head(mld)
  } ##- compute one MC draw and plot counts
  
  ##- apply MLD
  # if(!file.exists("./data/mc.mld.rds")){
    system.time(
      mc.mld <- replicate(Nsim, single.mc.mld(x = Winc))
    ) # 17s for Nsim = 100
  #  saveRDS(mc.mld, "./data/mc.mld.rds")
  # } else {
  #   mc.mld <- readRDS("../data/mc.mld.rds")
  # }
  ##- get rid of row = donor 'out'
  counts.mld <- mc.mld[-which(rownames(mc.mld)=='out'), ]
  return(list('states' = tip.states, 'counts' = counts.mld))
}

FILE = "./data/mc.mld.baseline.rds"
if(!file.exists(FILE)){
  system.time(
    list.counts.mld <- lapply(l, get.mc.mld.counts)
  ) # 1979s = 33mn for 100 sims * 100 MC reps
  saveRDS(list.counts.mld, FILE)
} else {
  list.counts.mld <- readRDS(FILE)
}

###################------

if (FALSE){ # stats
  PropOfNonZero <- lapply(list.counts.mld, function(x) round(apply(x$counts, 2, function(z) sum(z>0)) / dim(x$counts)[1] * 100, 2))
  s <- sapply(PropOfNonZero, summary) # 4% on non-zero
  m <- sapply(list.counts.mld, function(x) max(x$counts))
  table(m)
} ## stats counts

if(FALSE){
  ##---- 10 reps ----
  ##- test: bind 10 first MC replicates
  ks <- sample(1:100, 10)
  tab.mld.10 <- data.frame('id' = rep(rownames(counts.mld), length(ks)), 'mld' = as.vector(counts.mld[, ks]))
  table(tab.mld.10$mld)
  tab <- merge(x = tip.states, 
               y = tab.mld.10, by = 'id', all.x = T) # all.x = T)
  table(tab$mld, useNA = 'ifany')
  tab$mld[is.na(tab$mld)] <- 0L # all donors (if all.x = T)
  ##- test: bind 10 random MC replicates
  
  ##---- plots ----
  par(mfrow = c(1, 3))
  boxplot(stage ~ mld, data = tab, horizontal =TRUE, ylab = 'MLD count', xlab = 'stage')
  boxplot(age ~ mld, data = tab, horizontal = T, xlab= 'age sampling')
  boxplot(risk ~ mld, data=tab, horizontal = T, xlab = 'risk level')
} # plot on 10 reps

if(FALSE){
  ##---- first replicate ----
  ##- keep one df, merge copy MLD count
  k <- 3
  tab.mld <- data.frame('id' = rownames(counts.mld), 'mld' = counts.mld[, k])
  tab <- merge(x = tip.states,
               y = tab.mld, by = 'id', all.x = T) # all.y = T)
  tab$mld[is.na(tab$mld)] <- 0L # all donors (if all.x = T)
  
  ##---- models ----
  X4 <- c('factor(age)', 'factor(stage)', 'risk')
  Y = 'mld'
  Z = 1
  mod1 <- as.formula(paste(Y, '~', paste(X4, collapse = ' + ')))
  mod1bi <- as.formula(paste('binmld', '~', paste(X4, collapse = ' + ')))
  
  ##---- fit ----
  poisson1 <- glm(mod1 , family="poisson", data = tab) # poisson GLM
  nb1 <- glm.nb(mod1, data = tab) # negative binomial
  zi_poisson1 <- zeroinfl(mod1, data = tab, dist= 'poisson') # zero-inflated poisson
  zi_nb1 <- zeroinfl(mod1, data = tab, dist= 'negbin') # zero-inflated neg bin
  
  ## coefs
  fm <- list("PO" = poisson1, "NB" = nb1, "ZIPO" = zi_poisson1, "ZINB" = zi_nb1)
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
  
  tab.latex
  #html(tab.latex, rmarkdown = TRUE)
  
  ##- compute mean parameters across replicates and number of p-values < .05 by variables
  str(counts.mld)
  dimnames(counts.mld)
  str(tab.mld)
} ## tests zero-inflated regressions


## debug:  k = 1; counts = counts.mld; df = tip.states; X =  c('factor(age)', 'factor(stage)', 'risk')
get.mld.logit <- function(k, counts = counts.mld, df = tip.states, X =  c('factor(age)', 'factor(stage)', 'risk')){ # index column, states with ID variable
  tab.mld <- data.frame('id' = rownames(counts), 'mld' = counts[, k], stringsAsFactors = FALSE)
  tab <- merge(x = df,
               y = tab.mld, by = 'id', all.x = T) # all.y = T)
  tab$mld[is.na(tab$mld)] <- 0L # all donors (if all.x = T)
  
  ##- 0. logistic regression (0, non 0)
  mod1bi <- as.formula(paste('binmld', '~', paste(X, collapse = ' + ')))
  tab$binmld <- ifelse(tab$mld == 0, 0, 1)
  logit1 <- glm(mod1bi, family="binomial", data = tab)
  return(logit1)
}
# fits.logit <- lapply(1:dim(counts.mld)[2], function(z) get.mld.logit(k = z))
fits.logit <- lapply(list.counts.mld, function(x) lapply(1:dim(x$counts)[2], function(z) get.mld.logit(k = z, counts = x$counts, df = x$states)))

##- count any p-value < 0.05 by X
## debug: X = X; fit = fits.logit[[1]]
get.p.values <- function(fit, X =  c('factor(age)', 'factor(stage)', 'risk')){
  m <- summary(fit)
  p.values <- m$coef[,4] 
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
#p.sum <- sapply(fits.logit, get.p.values)
p.sum <- lapply(fits.logit, function(x) sapply(x, function(z) get.p.values(fit = z)))
p.sum.per.sim <- sapply(p.sum, rowSums)
apply(p.sum.per.sim, 1, summary)
