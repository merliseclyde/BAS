## ----setup, include=FALSE------------------------------------------------
#require(knitr)
require(MASS)

## ----data----------------------------------------------------------------
library(MASS)
data(UScrime)

## ----transform-----------------------------------------------------------
UScrime[,-2] = log(UScrime[,-2])

## ----bas-----------------------------------------------------------------
library(BAS)
crime.ZS =  bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null",
                   modelprior=uniform(), initprobs="eplogp") 

## ---- fig.show='hold'----------------------------------------------------
plot(crime.ZS, ask=F)


## ----pip, fig.width=5, fig.height=5--------------------------------------
plot(crime.ZS, which = 4, ask=FALSE, caption="", sub.caption="")

## ----print---------------------------------------------------------------
crime.ZS

## ----summary------------------------------------------------------------------
options(width = 80)
summary(crime.ZS)

## ----image, fig.width=5, fig.height=5-----------------------------------------
image(crime.ZS, rotate=F)

## ----coef---------------------------------------------------------------------
coef.ZS = coef(crime.ZS)


## ----plot---------------------------------------------------------------------
plot(coef.ZS, subset=c(5:6),  ask=F)

## ----coefall------------------------------------------------------------------

plot(coef.ZS, ask=FALSE)

## ----confint-coef-------------------------------------------------------------

confint(coef.ZS)

## ----plot-confint, fig.width=7------------------------------------------------
plot(confint(coef.ZS, parm=2:16))

## ---- warning=FALSE,  fig.width=7---------------------------------------------
plot(confint(coef(crime.ZS, estimator="HPM")))

## ---- warning=FALSE,  fig.width=7---------------------------------------------
plot(confint(coef(crime.ZS, estimator="MPM")))

## ----choice of estimator------------------------------------------------------
muhat.BMA = fitted(crime.ZS, estimator="BMA")
BMA  = predict(crime.ZS, estimator="BMA")

# predict has additional slots for fitted values under BMA, predictions under each model
names(BMA)

## ---- fig.width=5, fig.height=5-----------------------------------------------
par(mar=c(9, 9, 3, 3))
plot(muhat.BMA, BMA$fit, 
     pch=16, 
     xlab=expression(hat(mu[i])), ylab=expression(hat(Y[i])))
abline(0,1)

## ----HPM----------------------------------------------------------------------
HPM = predict(crime.ZS, estimator="HPM")

# show the indices of variables in the best model where 0 is the intercept
HPM$bestmodel

## -----------------------------------------------------------------------------
(crime.ZS$namesx[HPM$bestmodel +1])[-1]

## ----MPM----------------------------------------------------------------------
MPM = predict(crime.ZS, estimator="MPM")
(crime.ZS$namesx[attr(MPM$fit, 'model') +1])[-1]

## ----BPM----------------------------------------------------------------------
BPM = predict(crime.ZS, estimator="BPM")
(crime.ZS$namesx[attr(BPM$fit, 'model') +1])[-1]

## ---- fig.width=6, fig.height=6-----------------------------------------------
GGally::ggpairs(data.frame(HPM = as.vector(HPM$fit),  #this used predict so we need to extract fitted values
                   MPM = as.vector(MPM$fit),  # this used fitted
                   BPM = as.vector(BPM$fit),  # this used fitted
                   BMA = as.vector(BMA$fit))) # this used predict

## ----se, fig.width=7----------------------------------------------------------
BPM = predict(crime.ZS, estimator="BPM", se.fit=TRUE)
crime.conf.fit = confint(BPM, parm="mean")
crime.conf.pred = confint(BPM, parm="pred")
cbind(crime.conf.fit, crime.conf.pred)
plot(crime.conf.fit)

## ----pred---------------------------------------------------------------------

new.pred = predict(crime.ZS, newdata=UScrime, estimator="MPM")

## -----------------------------------------------------------------------------
system.time(
  for (i in 1:10) {
    crime.ZS <- bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null", method="BAS",
                   modelprior=uniform(), initprobs="eplogp") 
  }
)

system.time(
  for (i in 1:10)  {
    crime.ZS <-  bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null", method="deterministic",
                   modelprior=uniform(), initprobs="eplogp") 
  }
)


## ----MCMC---------------------------------------------------------------------
crime.ZS =  bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null",
                   modelprior=uniform(),
                   method = "MCMC") 

## ----diagnostics--------------------------------------------------------------
diagnostics(crime.ZS, type="pip",  pch=16)
diagnostics(crime.ZS, type="model",  pch=16)

## ----biggerMCMC, eval=FALSE---------------------------------------------------
#  crime.ZS =  bas.lm(y ~ .,
#                     data=UScrime,
#                     prior="ZS-null",
#                     modelprior=uniform(),
#                     method = "MCMC", MCMC.iterations = 10^6)
#  
#  # Don't run diagnostics(crime.ZS, type="model", pch=16)

## ----add-out------------------------------------------------------------------
data("stackloss")
stackloss$out.ind = diag(nrow(stackloss))

stack.bas = bas.lm(stack.loss ~ ., data=stackloss,
                method="MCMC", initprobs="marg-eplogp",
                prior="ZS-null", 
                modelprior=tr.poisson(4,10),
                MCMC.iterations=200000
               )

## -----------------------------------------------------------------------------
knitr::kable(as.data.frame(summary(stack.bas)))

