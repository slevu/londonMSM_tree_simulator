##---- libs2 ----
library(reshape2)
library("ggplot2")
theme_set(theme_bw())
library(scales)

##---- associations ----
## age vs sizes, degrees at different thr
# rm(list=ls())
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )

## change structure in long table (!!!! first few for speed)
.n <- 5 
.m <- sample(1:100, .n)
c <- lapply(cw_Baseline0, 
            function(x) do.call(rbind, x[.m]))

options(scipen = -100)
names(c)[-1] <- as.character( as.numeric(names(c)[-1]) )
options(scipen = 0) 
# rm(cw_Baseline0)

# head(c[[2]])
# names(c)
# str(c)
# lapply(c, function(x) mean(x$nbhsize))

##---- winsorize outliers -----
winsorize <- function (x, fraction=.05)
{
  if(length(fraction) != 1 || fraction < 0 ||
     fraction > 0.5) {
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction, 1-fraction), na.rm = TRUE)
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}

c_win <- lapply(c, function(x){
  if ( any(grepl('outdegree', names(x))) ) {
    x[, 'outdegree'] <- winsorize( x[, 'outdegree'] )
    return(x)
  } else {
    x[, 'size'] <- winsorize( x[, 'size'] )
    x[, 'nbhsize'] <- winsorize( x[, 'nbhsize'] )
    return(x)
  }
})

##---- reshape ----
molten <- function(a){
  .sa <- melt( a[[1]][, -c(6)],  measure.vars = 5)
# head(sa)
.thr <- melt(a[-1],  measure.vars = c(6, 7))
# head(thr)
tot <- rbind( cbind(.sa, "L1" = 'NA'), .thr[,-c(2,6)])
# head(tot)
return(tot)
}
tot <- molten(c)
tot_win <- molten(c_win)


##---- effect winsorize ----
tapply(tot$value, list(tot$L1, tot$variable), function(x) max(x, na.rm = T))
tapply(tot_win$value, list(tot_win$L1, tot_win$variable), function(x) max(x, na.rm = T))

##---- boxplot sizes vs age ---- 
p <- ggplot(tot, aes(x = factor(age), y = value, colour = variable), outlier.colour= alpha("black", 0.1)) + geom_boxplot() + theme(legend.position="none", legend.title = element_blank(), strip.text.x = element_text(size = 12)) 
# p2 <- p + facet_wrap(~ variable + L1, nrow = 4, scales = "free_y")
# p2
# p2 %+% tot_win
# p2 + geom_jitter(alpha = 0.1)
# 
# p3 <- ggplot(tot, aes(factor(age), value, fill = variable)) + geom_boxplot()
# p4 <- p3 + facet_wrap(~ variable + L1, nrow = 1, scales = "free_y")
# p4
# p4 %+% tot_win

##- at 0.015 only, winsorized
# p %+% tot[tot$L1 %in% c('NA','0.015'),] + facet_wrap(~ variable, nrow = 1, scales = "free_y")
p %+% tot_win[tot_win$L1 %in% c('NA','1.5e-02'),] + facet_wrap(~ variable + L1, nrow = 1, scales = "free_y", labeller=labeller(variable = c(outdegree = "Outdegree", size = "Cluster size", nbhsize = "Neighborhood \n size")))


##---- regressions ----
##- summarize regression on bootstrap
 
source('functions.R')
compare.reg.bs2

model1 <- "y ~ factor(age)"
model2 <- "y ~ factor(age) + factor(stage) + factor(age)*factor(stage)"

##---- m1 ----
m1 <- compare.reg.bs2(ls = cw_Baseline0, reg = lm, y = 'size', model = model1)
m1
##---- m2 ----
m2 <- compare.reg.bs2(ls = cw_Baseline0, reg = lm, y = 'size', model = model2)
m2
##---- m3 ----
m3 <- compare.reg.bs2(ls = cw_Baseline0, reg = lm, y = 'nbhsize', model = model2)
m3

####---- stop ----