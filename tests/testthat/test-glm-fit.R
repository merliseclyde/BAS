context("bayesglm fit")

# Github issue #67
test_that("bayesglm.fit", {
  data(Pima.tr, package="MASS")
  
  # Github issue #67

  expect_equal(bayesglm.fit(y = Pima.tr$type, 
                            x = Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = Pima.tr$type, 
                       x = Pima.tr$age, 
                       family=binomial())$coef)


  #OK below this
  
  expect_equal(bayesglm.fit(y = Pima.tr$type, 
                            x = Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = Pima.tr$type, 
                       x = Pima.tr$age, 
                       family=binomial())$coef)
  
  
  
  expect_equal(bayesglm.fit(y = cbind(Pima.tr$type, Pima.tr$type), 
                            x = cbind(1, 1.0*Pima.tr$age), 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind(Pima.tr$type, Pima.tr$type), 
                       x = cbind(1, 1.0*Pima.tr$age), 
                       family=binomial())$coef)  
  
  expect_equal(bayesglm.fit(y = as.double(Pima.tr$type == "Yes"), 
                            x = Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = as.double(Pima.tr$type == "Yes"), 
                       x = Pima.tr$age, 
                       family=binomial())$coef
  )
  
  expect_equal(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), 5.0), 
                            x = Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind((Pima.tr$type == "Yes"), 5.0), 
                            x = Pima.tr$age, 
                            family=binomial())$coef
  )
  expect_equal(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                            x = cbind(1, 1.0*Pima.tr$age), 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                       x = cbind(1, 1.0*Pima.tr$age), 
                       family=binomial())$coef)
  expect_equal(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                            x = Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                       x = Pima.tr$age, 
                       family=binomial())$coef)
  expect_equal(bayesglm.fit(y = cbind(Pima.tr$type, 2.0), 
                            x = cbind(1, 1.0*Pima.tr$age), 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind(Pima.tr$type, 2.0), 
                       x = cbind(1, 1.0*Pima.tr$age), 
                       family=binomial())$coef)
  
  wt = sample(1:10, size=nrow(Pima.tr), replace=TRUE)
  expect_equal(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                            x = cbind(1, 1.0*Pima.tr$age), weights = wt,
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                       x = cbind(1,1.0*Pima.tr$age), weights = wt,
                       family=binomial())$coef)
  expect_equal(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), wt), 
                            x = cbind(1, 1.0*Pima.tr$age), 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind((Pima.tr$type == "Yes"), wt), 
                       x = cbind(1,1.0*Pima.tr$age), 
                       family=binomial())$coef)
  expect_equal(bayesglm.fit(y = c((Pima.tr$type == "Yes"), rep(0.0, nrow(Pima.tr))), 
                            x = 1.0*c(Pima.tr$age, Pima.tr$age),
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = c((Pima.tr$type == "Yes"), rep(0.0, nrow(Pima.tr))), 
                       x = 1.0*c(Pima.tr$age, Pima.tr$age),
                       family=binomial())$coef
               )
  data(crabs, package = "glmbb")
  crabs.bas <- bas.glm(satell ~ color + spine + width + weight,
                           data = crabs,
                           family = poisson(),
                           betaprior = bic.prior(), modelprior = uniform(),
                           method = "BAS", 
                           n.models = 1, 
                           include.always = ~ color + spine + width + weight
                          )
  
  crabs.bayesglmfit = bayesglm.fit(crabs$satell, x = model.matrix( ~ color + spine + width + weight, data = crabs),
                                   family = poisson(), coefprior = bic.prior())
  crabs.glmfit = glm.fit(crabs$satell, x = model.matrix( ~ color + spine + width + weight, data = crabs),
                                   family = poisson())
  
  expect_equal(crabs.bayesglmfit$coef, as.numeric(crabs.glmfit$coef))
  expect_equal(crabs.bayesglmfit$coef, unlist(crabs.bas$mle))
  
 
  expect_no_error(bayesglm.fit(crabs$satell, x = model.matrix( ~ color + spine + width + weight, data = crabs),
                            family = poisson(), coefprior = bic.prior()))
  expect_error(bayesglm.fit(as.numeric(crabs$satell), x = model.matrix( ~ color + spine + width + weight, data = crabs),
                            family = poisson(), coefprior = EB.local()))
  
  
})

test_that("bayesglm.fit for Binomial Data", {
 
Gegevens <- data.frame(
  Jaar      = seq(from=2020,to=2024),
  Totaal = c(29,18,19,15,15),
  FeitenBeoordeling = c(14,12,13,7,7),
  Beoordeling = c(13,8,10,5,5),
  MentionFB = cbind(c(14,12,13,7,7), c(15,6,6,8,8)),
  MentionB  = cbind(c(13,8,10,5,5), c(16,10,9,10,10))
)

bin.mod <- glm(cbind(Gegevens$MentionFB.1,Gegevens$MentionFB.2) ~
                 Jaar, data=Gegevens,
               family=binomial(link = "logit"))
bin.fit = glm.fit(y = cbind(Gegevens$MentionFB.1,Gegevens$MentionFB.2),
                  x = cbind(1.0, Gegevens$Jaar), 
                  family=binomial(link = "logit"))
bayesglm.fit = bayesglm.fit(y = cbind(Gegevens$MentionFB.1,Gegevens$MentionFB.2),
                       x = cbind(1.0, Gegevens$Jaar), 
                       family=binomial(link = "logit"))
expect_equal(as.numeric(bin.mod$coefficients), bayesglm.fit$coefficients) #no names in bas.fit
})
