context("bayesglm fit")

# Github issue #67
test_that("bayesglm.fit", {
  data(Pima.tr, package="MASS")
  expect_error(bayesglm.fit(y = cbind((Pima.tr$type == "Yes"), 5.0), x = 1.0*Pima.tr$age, 
                            family=binomial(), coefprior=bic.prior(100)))
  
})
