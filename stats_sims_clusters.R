# rm(list=ls())

cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
cw_EqualStage0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.EqualStage0.rds" )

a <- cw_Baseline0[[2]]

## bind and distribution
# library(data.table)
# b <- rbindlist(a)

c <- lapply(cw_Baseline0, function(x){
  do.call("rbind", x)
})
names(c)
str(c)

hist(c[[2]]$size)
hist(c[["0.015"]]$outdegree)
hist(b$indegree)
sapply(cw_Baseline0, function(x){
  summary(sapply(x, function(x) sum(x$binclus) / length(x$binclus)))
})

##- merge cluster assignement and W_0 
names(cw_Baseline0[["0.015"]])
c <- cw_Baseline0[["0.015"]][[14]]

boxplot(c$size ~ c$stage)
boxplot(c$size ~ c$age)
boxplot(c$size ~ c$risk)
boxplot(c$outdegree ~ c$stage)
boxplot(c$outdegree ~ c$age)
boxplot(c$outdegree ~ c$risk)

t.test(size ~ risk, data = c)
t.test(outdegree ~ risk, data = c)

d <- cw_EqualStage0[["0.015"]][[14]]
t.test(size ~ risk, data = d)
t.test(outdegree ~ risk, data = d)

## https://rpubs.com/corey_sparks/27239


###--- stats cluster ---###
# listUKclus <- l_bs_uk ## without patients data

####---- n clusters 2 ----
## number of different clusters (counting size 1)
sapply(cw_Baseline0, function(x){
  summary(sapply(x, function(x) {
    length(unique(x$ClusterID) )
  }))
})

####---- membership ----
## proportion of cluster membership
sapply(cw_Baseline0, function(x){
  summary(sapply(x, function(x) sum(x$binclus) / length(x$binclus)))
})

####---- mean size 2 ----
## stats of mean size
sapply(cw_Baseline0, function(x){
  summary(sapply(x, function(x) mean(x$size)))
})

####---- median size 2 ----
## stats of median size
sapply(cw_Baseline0, function(x){
  summary(sapply(x, function(x) median(x$size)))
})

## stats of max size
sapply(cw_Baseline0, function(x){
  summary(sapply(x, function(x) max(x$size)))
})

lapply(cw_Baseline0, function(x) lapply(x, function(df) aggregate(df$size, by = list("stage" = df$stage), mean )))
head(cw_Baseline0[[1]][[1]])
####---- plots ----

tab_b <- l_Baseline0[["0.015"]][[40]]
tab_e <- l_EqualStage0[["0.015"]][[40]]
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
## type 2
mean(sapply(l_Baseline0[["0.015"]], u.test) > 0.05)
## type 1
mean(sapply(l_EqualStage0[["0.015"]], u.test) < 0.05)

## after downsampling

down_baseline <- lapply(l_Baseline0[["0.015"]], function(x) downsample(x, iter = 10))
down_equal <- lapply(l_EqualStage0[["0.015"]], function(x) downsample(x, iter = 10))

names(down_baseline[[1]])
mean(sapply(down_baseline, function(x) sapply(x, u.test)) > 0.05)
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