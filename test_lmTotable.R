### test output from lm in table ###

# require(papeR)
# library("devtools")
# install_github("hofnerb/papeR")

## ex. lm
## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
lm.D90 <- lm(weight ~ group - 1) # omitting intercept

fit <- summary(lm.D9)
a <- (capture.output(fit))

grep("Intercept", a)
grep("---", a)

b <- a[(grep("Intercept", a)+1):(grep("---", a)-1)]

DF <- read.table(textConnection( b ), fill = TRUE)
names(DF) <- c(" ", colnames(coef(summary(lm.D9))), " ")

# from https://stat.ethz.ch/pipermail/r-help/2008-July/166696.html
# DF <- read.table(textConnection(capture.output(summary(lm.D9))[11:12]), fill = TRUE)
# names(DF) <- c(" ", colnames(coef(summary(lm.D9))), " ")
# print(xtable(DF), include.rownames = FALSE)
