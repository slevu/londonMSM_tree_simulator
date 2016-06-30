###--- external code for stats_sims.Rmd
# rm(list=ls())

##---- libs ----
library(ggplot2)
library(cowplot)
library(reshape2)

##---- load data ----
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
cw <- cw_Baseline0[-c(2,3)]
names(cw)

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

##---- boxplot 1 ----
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
    facet_grid(~ method, scales = "free", space = "free", labeller=labeller  (method = c(SA = "SA", CL = "UCSD cluster", NB = "Neighborhood")))  + 
    xlab("Distance threshold") + ylab("p-value") + 
    theme(legend.position="none") +
    theme(strip.background = element_blank()) +
    background_grid()
  bp2
}

super_boxplot(p_uni)

##---- sign uni
##- table of proportion p-values < 0.05
tab1 <- data.frame("method" = names(p_uni),
                  "proportion univariate" = sapply(p_uni, function(x) mean( x < 0.05)), row.names = NULL)

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

super_boxplot(p_mult)

tab <- cbind(tab1, "proportion multivariate" = sapply(p_mult, function(x) mean( x < 0.05)))
row.names(tab) <- NULL
tab


###############------ surplus ------

##- essai neighborhood size
cor(l$size, l$nbhsize)
cor(l$size, l$risk)
cor(l$nbhsize, l$risk)

cor(l$outdegree, l$nbhsize, use = "complete" )
plot(l$outdegree, l$nbhsize)
boxplot(l$outdegree ~ l$risk) 
boxplot(l$nbhsize ~ l$risk) 
boxplot(l$size ~ l$risk) 

ggplot(l, aes(factor(risk), log(nbhsize))) + geom_violin()

###########################
###########################



##----  long ----
##- calculate individual mean size and outdegree over sims
## change structure in long table (first few for speed and 4 higher thr)
.n <- 3 
.m <- sample(1:100, .n)
a <- lapply(cw, function(x) do.call(rbind, x[.m]))
# str(a)

##- make sure id are ordered
sa <- list('risk' = as.character(a[[1]]$risk), 'SA' = a[[1]]$outdegree)
nb <- lapply(a[-1], function(x) c(x$nbhsize))
cl <- lapply(a[-1], function(x) c(x$size))
str(nb)
str(cl)
str(sa)

mu <- c(sa, 'CL' = cl, 'NB'= nb )
str(mu)
#####################-------------------------------laaaaaaa
head(mu)
mumu <- reshape(mu, varying = list(names(mu)[-1]), direction = 'long')
mumu <- melt(mu, id = 'risk')
str(mumu)
head(mumu)
?reshape
?melt
##- same order of id (already done)
# x <- lapply(a, function(df) df[order(df$id),] )

str(x)
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
    facet_grid(~ method, scales = "free", space = "free", labeller=labeller  (method = c(SA = "SA", CL = "UCSD cluster", NB = "Neighborhood")))  + 
    xlab("Distance threshold") + ylab("p-value") + 
    theme(legend.position="none") +
    theme(strip.background = element_blank()) +
    background_grid()
  bp2
}


##---- stats cluster ----
## summary
sapply(a[-1], function(x) summary(x$size))
## proportion in cluster
sapply(a[-1], function(x) mean(x$binclus))
## proportion in cluster by risk
sapply(a[-1], function(x) tapply(x$binclus, x$risk, mean))

##---- corr ----
##- correlation cluster size vs outdegree
for (i in 1:length(a)){
  print(paste(names(a)[i],
              round(cor(a[[i]]$size, 
                        a[[i]]$outdegree, 
                        use = "complete"), 3)) )
}
cor(a[["0.015"]]$size, a[["SA"]]$outdegree, use = "complete")
##---- plot long table ----
# quartz()
par(mfrow = c(length(a[-1])/2, 2), bty = 'n')
for (i in 1:length(a[-1])){
  hist(a[-1][[i]]$size, main = names(a[-1])[i], xlab = '')
}
dev.off()

##- hist log - log
par(mfrow = c(length(a)/2, 2), bty = 'n')
for (i in 1:length(a)){
  h <- hist(log(a[[i]]$size), plot = F)
  h$counts <- log1p(h$counts) # log(y)
  plot(h, ylab = "log(Freq)",
       main = names(a)[i], xlab = "log(size)")
}
dev.off()

##- size vs outdegree
par(mfrow = c(length(a)/2, 2), bty = 'n')
for (i in 1:length(a)){
  plot(a[[i]]$size, a[[i]]$outdegree,
       xlab = "cluster size",
       ylab = "out-degree",
       col="#00000050", 
       main = names(a)[i])
}
dev.off()

##---- combine into means ----
## mean by id of set of variables var
# mean.sims <- lapply(a, function(df) {
#   aggregate(df[, c("binclus", "size", "outdegree", "indegree")],
#             list("id" = df$id, 
#                  "stage" = df$stage,
#                  "age" = df$age,
#                  "risk" = df$risk),                              function(x) mean(x, na.rm = TRUE) )
# })



# 
# par(mfrow = c(1, 2), bty = 'n')
# hist(a[[1]]$outdegree, main = "", xlab = "out-degree")
# hist(a[[1]]$indegree, main = "", xlab =  "in-degree")
# 
# dev.off()
# 
# ##- mean
# str(mean.sims)
# b <- mean.sims
# 
# par(mfrow = c(length(b), 2), bty = 'n')
# for (i in 1:length(b)){
#   hist(b[[i]]$size, main = names(b)[i], xlab = "cluster size")
#   plot(b[[i]]$size, b[[i]]$outdegree,
#        xlab = "cluster size",
#        ylab = "out-degree",
#        col="#00000050", 
#        main = names(b)[i])
# }
# dev.off()

#####
#####
#####
##- reorganize: for each sim, merge SA and cluster sizes for different thresholds
revert_list <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list) 
}

##- revert and keep 4 thresholds
rcw <- revert_list(cw_Baseline0[-c(2,3)])
# names(cw); names(cw[[1]])

r <- rcw[[1]]

tabs <- lapply(cw, function(x) {
  cbind(x[["SA"]], x[["0.015"]][ match(x[["SA"]]$id, x[["0.015"]]$id), c('size','nbhsize') ])
})

str(tabs)

# p_cs2 <- lapply(tabs, function (x)  {
# u.test.risk(x, y = "size")
# } )