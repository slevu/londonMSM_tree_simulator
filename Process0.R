####---- run model0.R ----####
#### 
####---- process o and tree ----####
## Load newest Rdata
l <-  list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
load(l[length(l)])
str(o)
require(deSolve)
?ode
