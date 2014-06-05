require(MASS)
library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",
                    modelprior=beta.binomial(1,1),
                    initprobs= "eplogp") 
summary(crime.bic)
plot(crime.bic)
image(crime.bic, subset=-1)


# takes a while to run:
# crime.coef = coefficients(crime.bic)
# crime.coef
# par(mfrow=c(3,2))
# plot(crime.coef, ask=FALSE)

# see update
#crime.aic = update(crime.bic, newprior="AIC")
#crime.zs = update(crime.bic, newprior="ZS-null")

#crime.EBG = EB.global.bma(crime.bic)
# same as update(crime.bic, newprior="EB-global")
#image(crime.EBG, subset=-1)

