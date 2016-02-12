####---- run model0.R ----####
#### rm(list=ls())
####---- process o and tree ----####
## Load newest Rdata
l <-  list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
load(l[length(l)])
str(o)
require(deSolve)
?ode
##- ode(y, times, func, parms, method)
