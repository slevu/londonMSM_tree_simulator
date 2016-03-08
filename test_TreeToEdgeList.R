library(ape)
tt <- rcoal(n = 20)
# str(tt) plot(tt)
source("TreeToEdgeList.R")
TreeToEdgeList(tt, stats = F)

