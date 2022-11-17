context("bayesglm fit")

# Github issue #67
test_that("bayesglm.fit", {
  data(Pima.tr, package="MASS")
  
  # Github issue #67
  expect_error(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), 5.0), x = Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior()))
  
  expect_error(bayesglm.fit(y = Pima.tr$type, x = 1.0*Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior()))
  expect_error(bayesglm.fit(y = 1.0*(Pima.tr$type == "Yes"), x = Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior()))
  
  #OK below this
  expect_equal(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                            x = 1.0*Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = cbind((Pima.tr$type == "Yes"), 1.0), 
                       x = 1.0*Pima.tr$age, 
                       family=binomial())$coef)
  expect_equal(bayesglm.fit(y = c((Pima.tr$type == "Yes"), rep(0.0, nrow(Pima.tr))), 
                            x = 1.0*c(Pima.tr$age, Pima.tr$age),
                            family=binomial(), coefprior=bic.prior())$coef,
               glm.fit(y = c((Pima.tr$type == "Yes"), rep(0.0, nrow(Pima.tr))), 
                       x = 1.0*c(Pima.tr$age, Pima.tr$age),
                       family=binomial())$coef
               )
})
