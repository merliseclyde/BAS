loc = system.file("testdata", package="BAS")
d = read.csv(paste(loc, "JASP-testdata.csv", sep="/"))

simpleFormula = as.formula("contNormal ~ contGamma + contcor1 + contGamma * contcor1 ")

set.seed(1)
library(BAS)
basObj = bas.lm(simpleFormula,
                      data = d,
                      alpha = 0.125316,
                      prior = "JZS", include.always=as.formula("contNormal ~ contcor1"),
                      weights = d$facFifty)
traceback()
image(basObj, rotate=FALSE)
image(basObj, rotate=FALSE, drop.always.included=TRUE)
basObj$include.always
