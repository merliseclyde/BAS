pkgname <- "BAS"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BAS')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BAS")
### * BAS

flush(stderr()); flush(stdout())

### Name: BAS
### Title: BAS: Bayesian Model Averaging using Bayesian Adaptive Sampling
### Aliases: BAS BAS-package BAS-package
### Keywords: package regression

### ** Examples

data("Hald")
hald.gprior =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="g-prior")

# more complete demos
demo(BAS.USCrime)
demo(BAS.hald)




cleanEx()
nameEx("Bernoulli")
### * Bernoulli

flush(stderr()); flush(stdout())

### Name: Bernoulli
### Title: Independent Bernoulli Prior Distribution for Models
### Aliases: Bernoulli bernoulli

### ** Examples

Bernoulli(.9)




cleanEx()
nameEx("CCH")
### * CCH

flush(stderr()); flush(stdout())

### Name: CCH
### Title: Generalized g-Prior Distribution for Coefficients in BMA Models
### Aliases: CCH

### ** Examples

CCH(alpha=.5, beta=100, s=0) 




cleanEx()
nameEx("EB.global")
### * EB.global

flush(stderr()); flush(stdout())

### Name: EB.global
### Title: Find the global Empirical Bayes estimates for BMA
### Aliases: EB.global EB.global.bas
### Keywords: regression

### ** Examples


library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
# EB local uses a different g within each model
crime.EBL =  bas.lm(y ~ ., data=UScrime, n.models=2^15,
                    prior="EB-local", initprobs= "eplogp")
# use a common (global) estimate of g
crime.EBG = EB.global(crime.EBL)





cleanEx()
nameEx("EB.local")
### * EB.local

flush(stderr()); flush(stdout())

### Name: EB.local
### Title: Empirical Bayes Prior Distribution for Coefficients in BMA Model
### Aliases: EB.local EB

### ** Examples

EB.local()





cleanEx()
nameEx("IC.prior")
### * IC.prior

flush(stderr()); flush(stdout())

### Name: IC.prior
### Title: Information Criterion Families of Prior Distribution for
###   Coefficients in BMA Models
### Aliases: IC.prior aic.prior AIC.prior bic.prior BIC.prior

### ** Examples

IC.prior(2)
          aic.prior()
          bic.prior(100)
          



cleanEx()
nameEx("Jeffreys")
### * Jeffreys

flush(stderr()); flush(stdout())

### Name: Jeffreys
### Title: Jeffreys Prior Distribution for $g$ for Mixtures of g-Priors for
###   Coefficients in BMA Models
### Aliases: Jeffreys

### ** Examples

Jeffreys()





cleanEx()
nameEx("TG")
### * TG

flush(stderr()); flush(stdout())

### Name: TG
### Title: Generalized g-Prior Distribution for Coefficients in BMA Models
### Aliases: TG

### ** Examples


TG(alpha=2)
CCH(alpha=2, beta=100, s=0)




cleanEx()
nameEx("bas.glm")
### * bas.glm

flush(stderr()); flush(stdout())

### Name: bas.glm
### Title: Bayesian Adaptive Sampling Without Replacement for Variable
###   Selection in Generalized Linear Models
### Aliases: bas.glm
### Keywords: GLM regression

### ** Examples


library(MASS)
data(Pima.tr)

pima.cch = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7, method="BAS",
              betaprior=CCH(a=1, b=532/2, s=0), family=binomial(),
              modelprior=beta.binomial(1,1))

summary(pima.cch)
image(pima.cch)

pima.robust = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
              method="BAS",
              betaprior=robust(), family=binomial(),
              modelprior=beta.binomial(1,1))

pima.BIC = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
              method="BAS",
              betaprior=bic.prior(), family=binomial(),
              modelprior=uniform())





cleanEx()
nameEx("bas.lm")
### * bas.lm

flush(stderr()); flush(stdout())

### Name: bas.lm
### Title: Bayesian Adaptive Sampling Without Replacement for Variable
###   Selection in Linear Models
### Aliases: bas.lm bas
### Keywords: regression

### ** Examples


library(MASS)
data(UScrime)
crime.bic =  bas.lm(log(y) ~ log(M) + So + log(Ed) + 
                    log(Po1) + log(Po2) +
                    log(LF) + log(M.F) + log(Pop) + log(NW) +
                    log(U1) + log(U2) + log(GDP) + log(Ineq) +
                    log(Prob) + log(Time), 
                    data=UScrime, n.models=2^15, prior="BIC",
                    modelprior=beta.binomial(1,1),
                    initprobs= "eplogp") 
summary(crime.bic)
plot(crime.bic)
image(crime.bic, subset=-1)
# more complete demo's
demo(BAS.hald)
## Not run: demo(BAS.USCrime) 




cleanEx()
nameEx("bayesglm.fit")
### * bayesglm.fit

flush(stderr()); flush(stdout())

### Name: bayesglm.fit
### Title: Fitting Generalized Linear Models and Bayesian marginal
###   likelihood evaluation
### Aliases: bayesglm.fit
### Keywords: GLM regression

### ** Examples


require(MASS)
library(MASS)
Pima.tr
Y = as.numeric(Pima.tr$type) - 1
X = cbind(1, as.matrix(Pima.tr[,1:7]))
out = bayesglm.fit(X, Y, family=binomial(),coefprior=bic.prior(n=length(Y)))
out$coef
out$se
# using built in function
glm(type ~ ., family=binomial(), data=Pima.tr)





cleanEx()
nameEx("beta.binomial")
### * beta.binomial

flush(stderr()); flush(stdout())

### Name: beta.binomial
### Title: Beta-Binomial Prior Distribution for Models
### Aliases: beta.binomial Beta.Binomial

### ** Examples

beta.binomial(1,10) #' @family priors modelpriors




cleanEx()
nameEx("beta.prime")
### * beta.prime

flush(stderr()); flush(stdout())

### Name: beta.prime
### Title: Beta-Prime Prior Distribution for Coefficients in BMA Model
### Aliases: beta.prime

### ** Examples

beta.prime(n=100)




cleanEx()
nameEx("bodyfat")
### * bodyfat

flush(stderr()); flush(stdout())

### Name: bodyfat
### Title: Bodyfat Data
### Aliases: bodyfat Bodyfat
### Keywords: datasets

### ** Examples


data(bodyfat)
bodyfat.bas = bas.lm(Bodyfat ~ Abdomen, data=bodyfat, prior="ZS-null")
summary(bodyfat.bas)
plot(Bodyfat ~ Abdomen, data=bodyfat, xlab="abdomen circumference (cm)")
betas = coef(bodyfat.bas)$postmean   # current version has that intercept is ybar
betas[1] = betas[1] - betas[2]*bodyfat.bas$mean.x
abline(betas)
abline(coef(lm(Bodyfat ~ Abdomen, data=bodyfat)), col=2, lty=2)




cleanEx()
nameEx("coef")
### * coef

flush(stderr()); flush(stdout())

### Name: coef.bas
### Title: Coefficients of a Bayesian Model Average object
### Aliases: coef.bas coef coefficients coefficients.bas print.coef.bas
###   print.coef.bas
### Keywords: regression

### ** Examples


data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, n.models=2^4, alpha=13,
                      prior="ZS-null", initprobs="Uniform", update=10)
coef.hald.gprior = coefficients(hald.gprior)
coef.hald.gprior
plot(coef.hald.gprior)
confint(coef.hald.gprior)

#Estimation under Median Probability Model
coef.hald.gprior = coefficients(hald.gprior, estimator="MPM")
coef.hald.gprior
plot(coef.hald.gprior)
plot(confint(coef.hald.gprior))


coef.hald.gprior = coefficients(hald.gprior, estimator="HPM")
coef.hald.gprior
plot(coef.hald.gprior)
confint(coef.hald.gprior)

# To add estimation under Best Predictive Model





cleanEx()
nameEx("confint.coef")
### * confint.coef

flush(stderr()); flush(stdout())

### Name: confint.coef.bas
### Title: Compute Credible Intervals for BAS regression coefficients from
###   BAS objects
### Aliases: confint.coef.bas confint
### Keywords: regression

### ** Examples



data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
coef.hald = coef(hald.gprior)
confint(coef.hald)
confint(coef.hald, approx=FALSE, nsim=5000)





cleanEx()
nameEx("confint.pred")
### * confint.pred

flush(stderr()); flush(stdout())

### Name: confint.pred.bas
### Title: Compute Credible (Bayesian Confidence) Intervals for a BAS
###   predict object
### Aliases: confint.pred.bas
### Keywords: regression

### ** Examples


data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
hald.pred = predict(hald.gprior, estimator="BPM", predict=FALSE, se.fit=TRUE) 
confint(hald.pred, parm="mean") 
confint(hald.pred)  #default
hald.pred = predict(hald.gprior, estimator="BMA", predict=FALSE, se.fit=TRUE) 
confint(hald.pred)





cleanEx()
nameEx("cv.summary.bas")
### * cv.summary.bas

flush(stderr()); flush(stdout())

### Name: cv.summary.bas
### Title: Summaries for Out of Sample Prediction
### Aliases: cv.summary.bas
### Keywords: regression

### ** Examples



## Not run: 
##D library(foreign)
##D cognitive = read.dta("http://www.stat.columbia.edu/~gelman/arm/examples/child.iq/kidiq.dta")
##D cognitive$mom_work = as.numeric(cognitive$mom_work > 1)
##D cognitive$mom_hs =  as.numeric(cognitive$mom_hs > 0)
##D colnames(cognitive) = c("kid_score", "hs","iq", "work", "age")
##D 
##D set.seed(42)
##D n = nrow(cognitive)
##D test = sample(1:n, size=round(.20*n), replace=FALSE)
##D testdata =  cognitive[test,]
##D traindata = cognitive[-test,]
##D cog_train = bas.lm(kid_score ~ ., prior="BIC", modelprior=uniform(), data=traindata)
##D yhat = predict(cog_train, newdata=testdata, estimator="BMA", se=F)
##D cv.summary.bas(yhat$fit, testdata$kid_score)
## End(Not run)



cleanEx()
nameEx("diagnostics")
### * diagnostics

flush(stderr()); flush(stdout())

### Name: diagnostics
### Title: BAS MCMC diagnostic plot.
### Aliases: diagnostics

### ** Examples


library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.ZS =  bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null",
                   modelprior=uniform(),
                   method = "MCMC",
                   MCMC.iter = 1000)   # short run for the example
diagnostics(crime.ZS)





cleanEx()
nameEx("eplogprob")
### * eplogprob

flush(stderr()); flush(stdout())

### Name: eplogprob
### Title: eplogprob - Compute approximate marginal inclusion probabilities
###   from pvalues
### Aliases: eplogprob
### Keywords: regression

### ** Examples


library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
eplogprob(lm(y ~ ., data=UScrime))





cleanEx()
nameEx("eplogprob.marg")
### * eplogprob.marg

flush(stderr()); flush(stdout())

### Name: eplogprob.marg
### Title: eplogprob.marg - Compute approximate marginal inclusion
###   probabilities from pvalues
### Aliases: eplogprob.marg
### Keywords: regression

### ** Examples


library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
eplogprob(lm(y ~ ., data=UScrime))



cleanEx()
nameEx("fitted")
### * fitted

flush(stderr()); flush(stdout())

### Name: fitted.bas
### Title: Fitted values for a BAS BMA objects
### Aliases: fitted.bas fitted
### Keywords: regression

### ** Examples


data(Hald)
hald.gprior =  bas.lm(Y~ ., data=Hald, prior="ZS-null", initprobs="Uniform")
plot(Hald$Y, fitted(hald.gprior, estimator="HPM"))
plot(Hald$Y, fitted(hald.gprior, estimator="BMA", top=3))
plot(Hald$Y, fitted(hald.gprior, estimator="MPM"))
plot(Hald$Y, fitted(hald.gprior, estimator="BPM"))




cleanEx()
nameEx("force.heredity.bas")
### * force.heredity.bas

flush(stderr()); flush(stdout())

### Name: force.heredity.bas
### Title: Post processing function to force constraints on interaction
###   inclusion bas BMA objects
### Aliases: force.heredity.bas
### Keywords: regression

### ** Examples

data(Hald)
bas.hald = bas.lm(Y ~ .^2, data=Hald)
bas.hald.int = force.heredity.bas(bas.hald)
image(bas.hal.int)



cleanEx()
nameEx("g.prior")
### * g.prior

flush(stderr()); flush(stdout())

### Name: g.prior
### Title: Families of G-Prior Distribution for Coefficients in BMA Models
### Aliases: g.prior

### ** Examples

g.prior(100)




cleanEx()
nameEx("hyper.g")
### * hyper.g

flush(stderr()); flush(stdout())

### Name: hyper.g
### Title: Hyper-g-Prior Distribution for Coefficients in BMA Models
### Aliases: hyper.g

### ** Examples

hyper.g(alpha=.5) 





cleanEx()
nameEx("hyper.g.n")
### * hyper.g.n

flush(stderr()); flush(stdout())

### Name: hyper.g.n
### Title: Generalized hyper-g/n Prior Distribution for g for mixtures of
###   g-priors on Coefficients in BMA Models
### Aliases: hyper.g.n

### ** Examples

n = 500
hyper.g.n(alpha = 3, n=n)




cleanEx()
nameEx("hypergeometric1F1")
### * hypergeometric1F1

flush(stderr()); flush(stdout())

### Name: hypergeometric1F1
### Title: Confluent hypergeometric2F1 function
### Aliases: hypergeometric1F1
### Keywords: math

### ** Examples

hypergeometric1F1(11.14756, 0.5, 0.00175097)





cleanEx()
nameEx("hypergeometric2F1")
### * hypergeometric2F1

flush(stderr()); flush(stdout())

### Name: hypergeometric2F1
### Title: Gaussian hypergeometric2F1 function
### Aliases: hypergeometric2F1
### Keywords: math

### ** Examples

hypergeometric2F1(12,1,2,.65)




cleanEx()
nameEx("image")
### * image

flush(stderr()); flush(stdout())

### Name: image.bas
### Title: Images of models used in Bayesian model averaging
### Aliases: image.bas image
### Keywords: regression

### ** Examples


require(graphics)
data("Hald")
hald.ZSprior =  bas.lm(Y~ ., data=Hald,  prior="ZS-null")
image(hald.ZSprior, subset=-1)




cleanEx()
nameEx("intrinsic")
### * intrinsic

flush(stderr()); flush(stdout())

### Name: intrinsic
### Title: Intrinsic Prior Distribution for Coefficients in BMA Models
### Aliases: intrinsic

### ** Examples

n = 500;
 tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n)





cleanEx()
nameEx("list2matrix")
### * list2matrix

flush(stderr()); flush(stdout())

### Name: list2matrix.bas
### Title: Coerce a BAS list object into a matrix.
### Aliases: list2matrix.bas list2matrix
### Keywords: regression

### ** Examples


## Not run: 
##D library(MASS)
##D data(UScrime)
##D UScrime[,-2] = log(UScrime[,-2])
##D crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
##D                     initprobs= "eplogp") 
##D coef = list2matrix.bas(crime.bic, "ols")  # extract all ols coefficients
##D se = list2matrix.bas(crime.bic, "ols.se")
##D models = list2matrix.which(crime.bic)     #matrix of model indicators
##D models = which.matrix(crime.bic$which, crime.bic$n.vars)     #matrix of model indicators
## End(Not run)




cleanEx()
nameEx("list2matrix.which")
### * list2matrix.which

flush(stderr()); flush(stdout())

### Name: list2matrix.which
### Title: Coerce a BAS list object into a matrix.
### Aliases: list2matrix.which
### Keywords: regression

### ** Examples


## Not run: 
##D library(MASS)
##D data(UScrime)
##D UScrime[,-2] = log(UScrime[,-2])
##D crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
##D                     initprobs= "eplogp") 
##D coef = list2matrix.bas(crime.bic, "ols")  # extract all ols coefficients
##D se = list2matrix.bas(crime.bic, "ols.se")
##D models = list2matrix.which(crime.bic)     #matrix of model indicators
##D models = which.matrix(crime.bic$which, crime.bic$n.vars)     #matrix of model indicators
## End(Not run)




cleanEx()
nameEx("phi1")
### * phi1

flush(stderr()); flush(stdout())

### Name: phi1
### Title: Compound Confluent hypergeometric function of two variables
### Aliases: phi1
### Keywords: math

### ** Examples


# special cases
# Phi1(a, b, c, x=0, y) = 2F1(b, a; c, y)
phi1(1, 2, 1.5, 0, 1/100);
hypergeometric2F1(2, 1, 1.5, 1/100, log = FALSE)

# Phi1(a, b=0, c, x, y) = Phi(a, b, c, x, y=0) = 1F1(a, c, x) ## ??
phi1(1, 0, 1.5, 3, 1/100);
hypergeometric1F1(1, 1.5, 3, log = FALSE);




cleanEx()
nameEx("plot")
### * plot

flush(stderr()); flush(stdout())

### Name: plot.bas
### Title: Plot Diagnostics for an BAS Object
### Aliases: plot.bas
### Keywords: regression

### ** Examples


data(Hald)
hald.gprior =  bas.lm(Y~ ., data=Hald, prior="g-prior", alpha=13,
                      modelprior=beta.binomial(1,1),
                      initprobs="eplogp")

plot(hald.gprior)





cleanEx()
nameEx("plot.coef")
### * plot.coef

flush(stderr()); flush(stdout())

### Name: plot.coef.bas
### Title: Plots the posterior distributions of coefficients derived from
###   Bayesian model averaging
### Aliases: plot.coef.bas
### Keywords: regression

### ** Examples


## Not run: 
##D library(MASS)
##D data(UScrime)
##D UScrime[,-2] = log(UScrime[,-2])
##D crime.bic = bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC")
##D plot(coefficients(crime.bic), ask=TRUE)
## End(Not run)




cleanEx()
nameEx("plot.confint")
### * plot.confint

flush(stderr()); flush(stdout())

### Name: plot.confint.bas
### Title: Plot Bayesian Confidence Intervals
### Aliases: plot.confint.bas
### Keywords: bayesian regression

### ** Examples


data(Hald)
hald.ZS = bas.lm(Y ~ ., data=Hald, prior="ZS-null", modelprior=uniform())
plot(confint(coef(hald.ZS),parm=2:5))
plot(confint(predict(hald.ZS, se.fit=TRUE), parm="mean"))




cleanEx()
nameEx("predict")
### * predict

flush(stderr()); flush(stdout())

### Name: predict.bas
### Title: Prediction Method for an object of class BMA
### Aliases: predict.bas predict
### Keywords: regression

### ** Examples


data("Hald")
hald.gprior =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="g-prior")

predict(hald.gprior, newdata=Hald, estimator="BPM", se.fit=TRUE, prediction=FALSE)
# same as fitted
fitted(hald.gprior,estimator="BPM")

# default is BMA and estimation of mean vector
hald.bma = predict(hald.gprior, top=5, se.fit=TRUE)  
confint(hald.bma)

hald.BPM = predict(hald.gprior, newdata=Hald[1,],
                    prediction=TRUE, se.fit=TRUE,
                    estimator="BPM") 
confint(hald.BPM)

hald.hpm = predict(hald.gprior, newdata=Hald[1,],
                    prediction=TRUE, se.fit=TRUE,
                    estimator="HPM") 
confint(hald.hpm)

hald.mpm = predict(hald.gprior, newdata=Hald[1,],
                    prediction=TRUE, se.fit=TRUE,
                    estimator="MPM") 
confint(hald.mpm)




cleanEx()
nameEx("predict.basglm")
### * predict.basglm

flush(stderr()); flush(stdout())

### Name: predict.basglm
### Title: Prediction Method for an object of class basglm
### Aliases: predict.basglm
### Keywords: regression

### ** Examples


library(MASS)
data(Pima.tr)
data(Pima.te)
Pima.bas = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7, method="BAS",
           betaprior=CCH(a=1, b=nrow(Pima.tr)/2, s=0), family=binomial(),
           modelprior=uniform())
 pred = predict(Pima.bas, newdata=Pima.te, top=1)  # Highest Probability model
 cv.summary.bas(pred$fit, Pima.te$type, score="miss-class")




cleanEx()
nameEx("print")
### * print

flush(stderr()); flush(stdout())

### Name: print.bas
### Title: Print a Summary of Bayesian Model Averaging objects from BAS
### Aliases: print.bas print
### Keywords: print regression

### ** Examples


library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
print(crime.bic)
summary(crime.bic)



cleanEx()
nameEx("robust")
### * robust

flush(stderr()); flush(stdout())

### Name: robust
### Title: Robust-Prior Distribution for Coefficients in BMA Model
### Aliases: robust

### ** Examples

robust(100)





cleanEx()
nameEx("summary")
### * summary

flush(stderr()); flush(stdout())

### Name: summary.bas
### Title: Summaries of Bayesian Model Averaging objects from BAS
### Aliases: summary.bas summary
### Keywords: print regression

### ** Examples


library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
print(crime.bic)
summary(crime.bic)



cleanEx()
nameEx("tCCH")
### * tCCH

flush(stderr()); flush(stdout())

### Name: tCCH
### Title: Generalized tCCH g-Prior Distribution for Coefficients in BMA
###   Models
### Aliases: tCCH

### ** Examples

n = 500;
 tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n)



cleanEx()
nameEx("testBF.prior")
### * testBF.prior

flush(stderr()); flush(stdout())

### Name: testBF.prior
### Title: Test based Bayes Factors for BMA Models
### Aliases: testBF.prior

### ** Examples


testBF.prior(100)
library(MASS)
data(Pima.tr)

# use g = n 
bas.glm(type ~ ., data=Pima.tr, family=binomial(), 
        betaprior=testBF.prior(nrow(Pima.tr)),
        modelprior=uniform(), method="BAS")



cleanEx()
nameEx("tr.beta.binomial")
### * tr.beta.binomial

flush(stderr()); flush(stdout())

### Name: tr.beta.binomial
### Title: Truncated Beta-Binomial Prior Distribution for Models
### Aliases: tr.beta.binomial tr.Beta.Binomial

### ** Examples


tr.beta.binomial(1,10, 5)
library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
                    modelprior=tr.beta.binomial(1,1,8),
                    initprobs= "eplogp")





cleanEx()
nameEx("tr.poisson")
### * tr.poisson

flush(stderr()); flush(stdout())

### Name: tr.poisson
### Title: Truncated Poisson Prior Distribution for Models
### Aliases: tr.poisson tr.Poisson

### ** Examples

tr.poisson(10, 50)




cleanEx()
nameEx("tr.power.prior")
### * tr.power.prior

flush(stderr()); flush(stdout())

### Name: tr.power.prior
### Title: Truncated Power Prior Distribution for Models
### Aliases: tr.power.prior tr.Power.Prior

### ** Examples


tr.power.prior(2, 8)
library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
                    modelprior=tr.power.prior(2,8),
                    initprobs= "eplogp")





cleanEx()
nameEx("uniform")
### * uniform

flush(stderr()); flush(stdout())

### Name: uniform
### Title: Uniform Prior Distribution for Models
### Aliases: uniform Uniform

### ** Examples

uniform()




cleanEx()
nameEx("update")
### * update

flush(stderr()); flush(stdout())

### Name: update.bas
### Title: Update BAS object using a new prior
### Aliases: update.bas update
### Keywords: regression

### ** Examples


## Not run: 
##D library(MASS)
##D data(UScrime)
##D UScrime[,-2] = log(UScrime[,-2])
##D crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
##D crime.ebg = update(crime.bic, newprior="EB-global")
##D crime.zs = update(crime.bic, newprior="ZS-null")
## End(Not run)




cleanEx()
nameEx("which.matrix")
### * which.matrix

flush(stderr()); flush(stdout())

### Name: which.matrix
### Title: Coerce a BAS list object of models into a matrix.
### Aliases: which.matrix
### Keywords: regression

### ** Examples


## Not run: 
##D library(MASS)
##D data(UScrime)
##D UScrime[,-2] = log(UScrime[,-2])
##D crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
##D                     initprobs= "eplogp") 
##D models = which.matrix(crime.bic$which, crime.bic$n.vars)     #matrix of model indicators
## End(Not run)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
