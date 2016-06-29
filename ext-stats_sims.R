###--- external code for stats_sims.Rmd
# rm(list=ls())

##---- load data ----
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
# cw_EqualStage0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.EqualStage0.rds" )

##- reorganize: for each sim, merge SA and cluster sizes for different thresholds
revert_list <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list) 
}

cw <- revert_list(cw_Baseline0[-c(2,3)])
# names(cw); names(cw[[1]])
# ########################---------------- laaaaaa
r <- cw[[1]]
str(r)
system.time(a <- merge(r[["SA"]], r[["0.015"]], by = 'id') )
system.time(b <- cbind(r[["SA"]], r[["0.015"]][ match(r[["SA"]]$id, r[["0.015"]]$id), c('size','nbhsize') ]))
str(a)
str(b)
str(x)
test <- lapply(r[-1], function(x) cbind(x[ , c('size','nbhsize')] ))
test <- unlist(r[-1])
str(test)

test <- lapply(cw, function(x) {
  merge(x[["SA"]], x[["0.015"]], by = 'id')
})

system.time(
    z <- lapply(l_Baseline0, function(x){
      sapply(x, function(m) {
      # faster than `merge(df, s, all.x = TRUE)`
      #- add covariates
       # m <-  cbind(df, s[ match(df$id, s$id), -1 ])
      #- mean y by x
        agg <-  tapply(m$size, m$stage, mean )
        return(agg)
        }) 
      })
    )
  
  
##----  long ----
##- calculate individual mean size and outdegree over sims
## change structure in long table (first few for speed and 4 higher thr)
.n <- 3 
.m <- sample(1:100, .n)
a <- lapply(cw_Baseline0[-c(2,3)], 
             function(x) do.call(rbind, x[.m]))
# str(a)

##---- stats cluster ----
## summary
sapply(a[-1], function(x) summary(x$size))
## proportion in cluster
sapply(a[-1], function(x) mean(x$binclus))

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
par(mfrow = c(length(a)/2, 2), bty = 'n')
for (i in 1:length(a)){
  hist(a[[i]]$size, main = names(a)[i], xlab = '')
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

##---- u test ----
##- U tests on size or outdegree
u.test.risk <- function(df, y){
  U <- wilcox.test(df[df$risk == 2, y], 
                   df[df$risk == 1, y], 
                   alternative = "greater") # "two.sided", "less"
  return(U$p.value)
 }

##- p-values for test of cluster size in risk 2 vs risk 1
p_cs <- lapply(cw_Baseline0, function (x)  {
  sapply(x, function(df) u.test.risk(df, y = "size") )
} )
##- Add p-values for outdegree testing
p_sa <- sapply(cw_Baseline0[[1]], function(df) u.test.risk(df, y = "outdegree") )
p_uni <- c(p_cs, "SA"= list(p_sa) )

##- stats of p-values
sapply(p_uni, summary)

##- table of proportion p-values < 0.05
tab <- data.frame("method" = c(paste("Cluster", names(cw_Baseline0)), "SA"),
                  "proportion univariate" = sapply(p_uni, function(x) mean( x < 0.05)),
                  row.names = NULL)


##---- lm ----
##- lm model adjusting for stage of infection
# y = "outdegree"
lm.risk <- function(df, y){
  model0 <- "scale(log(y)) ~ scale(risk) * as.factor(stage)"
  model <- sub("y", y, model0)
  ##- pvalue for risk parameter
  p <- coef(summary( lm( model , data = df) ))[2, 4]
  return(p)
}
p_mult_cs <- lapply(cw_Baseline0, function (x)  {
  sapply(x, function(df) lm.risk(df, y = "size") )
} )
p_mult_sa <- sapply(cw_Baseline0[[1]], function(df) lm.risk(df, y = "outdegree") )
p_mult <- c(p_mult_cs, "SA"= list(p_mult_sa) )
sapply(p_mult, summary)

tab <- cbind(tab, "proportion multivariate" = sapply(p_mult, function(x) mean( x < 0.05)))
row.names(tab) <- NULL


##---- plot p-values
cc <- matrix(sapply(p_uni, c), ncol = length(p_uni))
colnames(cc) <- names(p_uni)
# head(cc)
d <- as.data.frame.table(cc)
# head(d)
# str(d)

ccc <- matrix(sapply(p_mult, c), ncol = length(p_mult))
colnames(ccc) <- names(p_mult)
# head(ccc)
dd <- as.data.frame.table(ccc)
# head(dd)
# str(dd)
library(ggplot2) 
##- univariate
ggplot(d[d$Var2 != "SA",], aes(x = Var2, y = Freq, fill = Var2)) + geom_violin()
##- multivariate
g <- ggplot(dd, aes(x = Var2, y = Freq, fill = Var2)) + geom_violin()
g

##---- compute age matrix for cluster size  ----
l <- cw_Baseline0[[2]][[2]]
str(l)
head(l[order(l$ClusterID),])

x <- l$ClusterID
names(x) <- l$age
source("functions.R")
el <- EdgeList(x)
head(el)
dim(el)
agemat <- matrix( 0, nrow = length(unique(el[,'from'])), ncol = length(unique(el[,'to'])) )
for (k in 1:dim(el)[1]){
    from <- el[k, 1]
    to <- el[k, 2]
    agemat[from, to] <- agemat[from, to] + 1
}
sum(agemat)
require(lattice)
# lattice.options(default.theme = standard.theme(color=F))
levelplot( mat2assortmat( agemat ), col.regions = heat.colors)

##- essai neighborhood size
cor(l$size, l$nbhsize)
cor(l$outdegree, l$nbhsize, use = "complete" )
plot(l$outdegree, l$nbhsize)
boxplot(l$outdegree ~ l$risk) 
boxplot(l$nbhsize ~ l$risk) 
boxplot(l$size ~ l$risk) 
