##---- import ----
## from Rcolgem
rm(list=ls())
w0 <- read.csv("W0.sim.csv")
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
y <- readRDS("demo.rds")
head(y)

##---- merge ----
out <-   merge(out0, y, 
               by.x = "patient", 
               by.y = "patient", 
               all.x = T, sort = FALSE)
# saveRDS(out, file = "../data/outdegree.rds")
head(out)

##---- plot1 ----
##- boxplot
boxplot(outdegree ~ age, data = out, xlab = "age",
        ylab = "outdegree")
boxplot(outdegree ~ stage, data = out, xlab = "stage",
        ylab = "outdegree")
boxplot(outdegree ~ risk, data = out, xlab = "risk",
        ylab = "outdegree")

##---- plot2 ----
##- continuous
plot(out$time / 365, out$outdegree)
     
##---- lm ----
##- continuous
summary(lm(outdegree ~ factor(age), data = out))
summary(lm(outdegree ~ factor(stage), data = out))
summary(lm(outdegree ~ risk, data = out))
summary(lm(outdegree ~ time, data = out))


##---- multivariate ----
##- continuous dob
# names(out_last3)
summary(lm(scale(outdegree) ~ scale(age) + scale(stage) + scale(risk) + scale(time) , data = out)) 
summary(lm(outdegree ~ factor(age) + factor(stage) + factor(risk) + time , data = out)) 





