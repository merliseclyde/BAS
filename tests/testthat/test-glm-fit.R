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
})
