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
