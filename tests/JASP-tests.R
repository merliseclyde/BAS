loc = system.file("testdata", package="BAS")
d = read.csv(paste(loc, "JASP-testdata.csv", sep="/"))

simpleFormula = as.formula("contNormal ~ contGamma + contcor1 + contGamma * contcor1 ")

library(BAS)
set.seed(1)
basObj = bas.lm(simpleFormula,
                data = d,
                alpha = 0.125316,
                prior = "JZS",
                include.always=as.formula("contNormal ~ contcor1"),
                modelprior=beta.binomial(1,1),
                weights = d$facFifty)

image(basObj, rotate=FALSE)
image(basObj, rotate=FALSE, drop.always.included=TRUE)
basObj$include.always

## old
##
## install.packages("BAS")

## library(BAS)

 set.seed(1)
 basObj.old = bas.lm(simpleFormula,
                data = d,
                alpha = 0.125316,
                prior = "JZS",
                include.always=as.formula("contNormal ~ contcor1"),
                modelprior=beta.binomial(),
                weights = d$facFifty, force.heredity = FALSE)
 basObj.old = force.heredity.bas(basObj.old)

 basObj.old$postprobs  #(check order of models)
 basOBJ$postprobs
