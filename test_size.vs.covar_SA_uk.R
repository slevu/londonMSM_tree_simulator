#### rm(list=ls())
#### 
source("functions.R")
       
out <-  readRDS(file = "../phylo-uk/data/outdegree.rds")

fit <- lm(outdegree ~ agediag + sqrt(cd4) + ydiag , data = out) 
summary(fit)
str(out)
### plot (input l as a list)
size.vs.covar(l = list(out), depvar = "outdegree", indepvar = c("agediag", "cd4", "ydiag"))

### reg
lm_model <- "outdegree ~ agediag + sqrt(cd4) + ydiag"
reg.sum(ls = list(out), reg = lm, model = lm_model)

