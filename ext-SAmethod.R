##---- import ----
## from Rcolgem
rm(list=ls())
w0 <- read.csv("../data/W0.csv")
head(w0[order(w0$donor),])
head(w0[order(w0$recipient),])

##---- symetric ----
## W_ij matrix is symetric ? yes
identical(w0$donor[order(w0$donor)], 
          w0$recipient[order(w0$recipient)])

##---- outdegree ----
## by patient
out0 <- aggregate(x = list(outdegree = w0$W),
                  by = list(patient = w0$donor), FUN = sum)
in0 <- aggregate(x = list(indegree = w0$W),
                  by = list(patient = w0$recipient), FUN = sum)
head(out0)
summary(out0$outdegree)
summary(in0$indegree)
hist(out0$outdegree)

##---- explanatory ----
## add individual explanatory variates
load("../data/sub.RData")
#head(df)
#names(df)
##- selection of df covariates
y <- df[,c("seqindex","patientindex", "dob_y",
           "agediag", "cd4", "vl", "onartflag",
           "ydiag", "agediag_cut", "cd4cut",
           "ydiag_cut", "CHICflag", "status")]
# str(y)

##---- merge ----
out <-   merge(out0, y, 
               by.x = "patient", 
               by.y = "seqindex", 
               all.x = T, sort = FALSE)
# saveRDS(out, file = "../data/outdegree.rds")

##---- quantiles ----
##- cut in 10 quantiles of dob
out$dob_qt <-  cut( out$dob_y, 
                  quantile(out$dob_y, prob = seq(0, 1, length = 11)) )

##---- last3 ----
## restricting to last 3 years
#table(out$ydiag)
out_last3 <-out[out$ydiag > 2007,]
hist(out_last3$outdegree)

##---- plot1 ----
##- quantiles
plot(out_last3$dob_qt, out_last3$outdegree)

##---- plot2 ----
##- continuous
plot(out_last3$dob_y, out_last3$outdegree)
     
##---- lm ----
##- continuous dob
lm_dob <- lm(outdegree ~ dob_y, data = out_last3) 
summary(lm_dob)

##- quantiles dob
lm_dob_qt <- lm(outdegree ~ dob_qt, data = out_last3) 
summary(lm_dob)

##---- multivariate ----
##- continuous dob
# names(out_last3)
lm <- lm(outdegree ~ agediag + sqrt(cd4) + ydiag , data = out_last3) 
summary(lm)

##- standardized
lms <- lm(scale(outdegree) ~ 
                 scale(agediag) + scale(sqrt(cd4)) + scale(ydiag) ,
                 data = out_last3) 
summary(lms)



##---- split periods of 3 years for univariate analyses of out-degree ----



