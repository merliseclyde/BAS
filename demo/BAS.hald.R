data(Hald)
hald.gprior =  bas.lm(Y~ ., data=Hald, prior="g-prior", alpha=13,
                      modelprior=beta.binomial(1,1),
                      initprobs="eplogp")

hald.gprior
plot(hald.gprior)
summary(hald.gprior)
image(hald.gprior, subset=-1, vlas=0)

hald.coef = coefficients(hald.gprior)
hald.coef
plot(hald.coef)
predict(hald.gprior, top=5, se.fit=TRUE)  
confint(predict(hald.gprior, Hald, estimator="BMA", se.fit=TRUE, top=5), parm="mean")
predict(hald.gprior, estimator="MPM", se.fit=TRUE)  
confint(predict(hald.gprior, Hald, estimator="MPM", se.fit=TRUE), parm="mean")

fitted(hald.gprior, estimator="HPM")
hald.gprior =  bas.lm(Y~ ., data=Hald, n.models=2^4,
                      prior="g-prior", alpha=13, modelprior=uniform(),
                      initprobs="eplogp")
hald.EB = update(hald.gprior, newprior="EB-global")
hald.bic = update(hald.gprior,newprior="BIC")
hald.zs = update(hald.bic, newprior="ZS-null")
