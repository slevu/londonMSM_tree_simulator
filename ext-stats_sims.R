###--- external code for stats_sims.Rmd
# rm(list=ls())

##---- libs ----
library(ggplot2)
library(cowplot)
library(reshape2)

##---- load data ----
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
cw <- cw_Baseline0[c('SA', '0.001', '0.005', '0.015', '0.05')]
# names(cw)

##- same order of id
cw <- lapply(cw, function(x){
  lapply(x, function(df) df[order(df$id),] )
})

##---- u test ----
##- U tests on size or outdegree
u.test.risk <- function(df, y){
  U <- wilcox.test(df[df$risk == 2, y], 
                   df[df$risk == 1, y], 
                   alternative = "greater") # "two.sided", "less"
  return(U$p.value)
 }

##- p-values for test of cluster size in risk 2 vs risk 1
p_cs <- lapply(cw[-1], function (x)  {
  sapply(x, function(df) u.test.risk(df, y = "size") )
} )
p_nb <- lapply(cw[-1], function (x)  {
  sapply(x, function(df) u.test.risk(df, y = "nbhsize") )
} )

##- p-values for outdegree testing
p_sa <- sapply(cw[[1]], function(df) u.test.risk(df, y = "outdegree") )

##- stats of p-values
p_uni <- c('CL' = p_cs, 'NB' = p_nb, "SA"= list(p_sa) )
# str(p_uni)

##- table of proportion p-values < 0.05
tab1 <- data.frame("method" = names(p_uni),
                  "univariate" = sapply(p_uni, function(x) mean( x < 0.05)), row.names = NULL)

##---- lm ----
##- lm model adjusting for stage of infection with option transformation
# y = "size"; df = cw[[4]][[2]]
lm.risk <- function(df, y){
  model0 <- "scale(y) ~ scale(risk) * as.factor(stage)"
  model <- sub("y", y, model0)
  ##- pvalue for risk parameter
  p <- coef(summary( lm( model , data = df) ))[2, 4]
  return(p)
}

p_mult_cs <- lapply(cw[-1], function (x)  {
  sapply(x, function(df) lm.risk(df, y = "size") )
} )
p_mult_nb <- lapply(cw[-1], function (x)  {
  sapply(x, function(df) lm.risk(df, y = "nbhsize") )
} )
p_mult_sa <- sapply(cw[[1]], function(df) {
  lm.risk(df[df$outdegree > 0,], y = "log(outdegree)") })

p_mult <- c('CL' = p_mult_cs, 'NB'= p_mult_nb, "SA"= list(p_mult_sa) )

tab <- cbind(tab1, "multivariate" = sapply(p_mult, function(x) mean( x < 0.05)))
row.names(tab) <- NULL

##---- logistic regression ----
# y = "size"; df = cw[[2]][[2]]
logit.risk <- function(df){
  model0 <- "binclus ~ factor(risk) + factor(age) + factor(stage)"
  ##- pvalue for risk parameter
  fit <- glm( model0 , data = df, family = binomial(link = 'logit'))
  cf <- coef( summary( fit ) )
  p <- cf[2, 4]
  effect <- cf[2,1]
  # or <- exp(cbind(OR = coef(fit), confint(fit)))[2, ] # too long with CI
  or <- unname(exp(coef(fit))[2])
  return(c(p = p, OR = or))
}
##- save p values and OR + CI (takes time !)
if(!file.exists("data/sim_ucsd_results2/logit_cs.rds")){
  system.time(
  logit_cs <- lapply(cw[-1], function (x)  {
    sapply(x, function(df) logit.risk(df) )
  } )
  ) # 34s
  saveRDS(logit_cs, file = "data/sim_ucsd_results2/logit_cs.rds")
} else {
  logit_cs <- readRDS(file = "data/sim_ucsd_results2/logit_cs.rds")
}
# str(logit_cs[[1]])
##- table of proportion p-values < 0.05 and mean OR
tab_logit <- data.frame("method" = names(logit_cs),
                   "p signif" = sapply(logit_cs, function(x) mean( x[1] < 0.05)),
                   "OR" = sapply(logit_cs, function(x) mean( x[2])) ,
                   row.names = NULL)


##---- boxplot risk 1 ----
super_boxplot <- function(ls){
  ##- long table
  a <- melt(ls)
  ## extra columns for facets
  a$method <- substr(a$L1, 1, 2)
  ##- add 'NA' for threshold of SA method
  a$thr <- c(regmatches(a$L1, regexpr("\\d\\..*|\\d*e[-+]?.*",  a$L1)), rep('NA', length(a[a$method == 'SA', 'method'])) )
  
  ##- plot
  bp1 <- ggplot(a, aes(reorder(thr, as.numeric(thr)), value, color = method)) + geom_boxplot()
  bp2 <- bp1  +
    facet_grid(~ method, scales = "free", space = "free", labeller=labeller  (method = c(SA = "SA", CL = "Cluster", NB = "Neighborhood")))  + 
    xlab("Distance threshold") + ylab("p-value") + 
    theme(legend.position="none") +
    theme(strip.background = element_blank()) +
    background_grid() +
    geom_hline(aes(yintercept = 0.05), linetype = "dashed")
  bp2
}

super_boxplot(ls = p_uni)

##---- boxplot risk 2 ----
super_boxplot(p_mult)

##---- boxplot base risk 1 ----
super_boxplot_base <- function(ls){
  # ls = p_uni
  ##- long table
  a <- melt(ls)
  ## extra columns for facets
  a$method <- substr(a$L1, 1, 2)
  ##- add 'NA' for threshold of SA method
  a$thr <- c(regmatches(a$L1, regexpr("\\d\\..*|\\d*e[-+]?.*",  a$L1)), rep('NA', length(a[a$method == 'SA', 'method'])) )
  # str(a); table(a$thr)
  
  ##- plot
  ##- to keep same range
  lim <- range(a$value)
  alpha <- 0.05
  # par(mfrow=c(1,3) ) #, bty = 'n')   
  layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), 
         widths=c(1,1.5,1.5), heights=c(1,1,1))
  
  boxplot(a[a$method == 'SA',]$value ~ a[a$method == 'SA',]$thr, ylim = lim)
  title(main = 'SA', font.main = 1, ylab = 'p-value', cex.lab = 1.2)
  abline(h = alpha, lty = 3)
  
  boxplot(a[a$method == 'CL',]$value ~ a[a$method == 'CL',]$thr, ylim = lim, yaxt="n")
  title(main = 'Cluster', xlab = 'Distance threshold', cex.lab = 1.2,  font.main = 1)
  abline(h = alpha, lty = 3)
  
  boxplot(a[a$method == 'NB',]$value ~ a[a$method == 'NB',]$thr, ylim = lim, yaxt="n")
  title(main = 'Neighborhood', xlab = 'Distance threshold', cex.lab = 1.2,  font.main = 1)
  abline(h = alpha, lty = 3)
}

super_boxplot_base(ls = p_uni)

##---- boxplot base risk 2 ----
super_boxplot_base(p_mult)


##----  long ----
##- calculate individual mean size and outdegree over sims
## change structure in long table (first few for speed and 4 higher thr)
.n <- 10 
.m <- sample(1:100, .n)
z <- lapply(cw, function(x) do.call(rbind, x[.m]))
# str(a)

##- make sure id are ordered
sa <- list('risk' = as.character(z[[1]]$risk), 'SA' = z[[1]]$outdegree)
nb <- lapply(z[-1], function(x) c(x$nbhsize))
cl <- lapply(z[-1], function(x) c(x$size))
# str(nb); str(cl); str(sa)

mu <- data.frame(sa, 'CL' = cl, 'NB'= nb )
# str(mu)

a <- melt(mu, id.vars = 'risk')

##---- boxplot 3 ----
# str(a); head(a); tail(a)
a$method <- substr(a$variable, 1, 2)

##- censored 5% outliers
source('functions.R')
a$value_win <- unlist(tapply(a$value, a$variable, function(x) (winsorize(x, 0.05)) ) )
# str(a); tail(a, 100)

##- plot
bp4  <- ggplot(a, aes(risk, value_win, color = method)) + geom_boxplot()
bp5 <- bp4  +
  facet_wrap(~ variable, scales = "free")  + 
  xlab("Risk level") + ylab("Outdegree or cluster size (censored)") + 
  theme(legend.position="none") +
  theme(strip.background = element_blank()) +
  background_grid()
bp5
# bp5 %+% a[a$value_win > 0,] 

##---- boxplot 4 ----
##- plot log, uncensored
bp6  <- ggplot(a, aes(risk, log(value + 10-5), color = method)) + geom_boxplot()
bp7 <- bp6  +
  facet_wrap(~ variable, scales = "free")  + 
  xlab("Risk level") + ylab("Outdegree or cluster size (log)") + 
  theme(legend.position="none") +
  theme(strip.background = element_blank()) +
  background_grid()
bp7

##---- stats cluster 1 ----
##- proportion in cluster by risk
round(sapply(z[-1], function(x) tapply(x$binclus, x$risk, mean)), 2)

##---- stats cluster 2 ----
##- mean outdegree or size, by risk
mm <- tapply(a$value, list(a$risk, a$variable), mean)
round(mm, 2)

##---- stats cluster 3 ----
##- ratio mean
round( mm[2,] / mm[1,], 2)

##------ stop ------
