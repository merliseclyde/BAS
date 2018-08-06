context("bas.lm")

test_that("shrinkage is less than or equal to 1", {
  data(Hald)
  hald.bas = bas.lm(Y ~ ., prior="ZS-null", modelprior=uniform(), data=Hald)
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas = bas.lm(Y ~ ., prior="EB-local", modelprior=uniform(), data=Hald)
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas = bas.lm(Y ~ ., prior="EB-global", modelprior=uniform(), data=Hald)
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas = bas.lm(Y ~ ., prior="hyper-g", modelprior=uniform(), data=Hald)
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas = bas.lm(Y ~ ., prior="hyper-g-n", modelprior=uniform(), data=Hald)
  expect_equal(0, sum(hald.bas$shrinkage > 1))

})

test_that("shrinkage is less than or equal to 1", {
  data(Hald)
  #  hald.bas = bas.lm(Y ~ ., prior="JZS", modelprior=uniform(), data=Hald)
  #  expect_equal(0, sum(hald.bas$shrinkage > 1))
})

test_that("A/BIC: shrinkage is equal to 1", {
  data(Hald)
  hald.BIC = bas.lm(Y ~ ., prior="BIC", modelprior=uniform(), data=Hald)
  expect_equal(hald.BIC$n.model, sum(hald.BIC$shrinkage == 1))
  hald.AIC = bas.lm(Y ~ ., prior="AIC", modelprior=uniform(), data=Hald)
  expect_equal(hald.AIC$n.model, sum(hald.AIC$shrinkage == 1))

})

test_that("no method", {
  data(Hald)
  expect_error(bas.lm(Y ~ ., prior="garbage",
                      modelprior=uniform(), data=Hald))
})

test_that("deterministic, BAS and MCMC+BAS", {
  data(Hald)
  hald.bas = bas.lm(Y ~ ., prior="BIC",
                      modelprior=uniform(), data=Hald)
  hald.MCMCbas = bas.lm(Y ~ ., prior="BIC", method="MCMC+BAS",
                    modelprior=uniform(), data=Hald, MCMC.iterations=1000)
  hald.deterministic = bas.lm(Y ~ ., prior="BIC",
                              method="deterministic",
                              modelprior=uniform(), data=Hald)
  expect_equal(hald.bas$probne0, hald.deterministic$probne0)
  expect_equal(hald.bas$probne0, hald.MCMCbas$probne0)
})

test_that("pivot", {
  data(Hald)
  hald.bas = bas.lm(Y ~ ., prior="BIC",
                    modelprior=uniform(), data=Hald)
  hald.deterministic = bas.lm(Y ~ ., prior="BIC",
                              method="deterministic",
                              modelprior=uniform(), data=Hald, pivot=TRUE)
  expect_equal(hald.bas$probne0, hald.deterministic$probne0)
})


test_that("pivoting with non-full rank design", {
  set.seed(42)
  dat = data.frame(Y = rnorm(5), X1=1:5, X2=1:5, X3 = rnorm(5))

  tmp.bas = bas.lm(Y ~ ., data=dat, prior="BIC", modelprior=uniform(), method="BAS", pivot=T)

  tmp.mcmc = bas.lm(Y ~ ., data=dat, prior="BIC", modelprior=uniform(), method="MCMC", pivot=T, MCMC.iterations=10000)
  expect_equal(sort(tmp.bas$R2), sort(tmp.mcmc$R2))
})

test_that("prediction versus fitted", {
  data(Hald)
  hald.ZS = bas.lm(Y ~ ., prior="ZS-null", modelprior=uniform(), data=Hald)
  expect_equal(as.vector(fitted(hald.ZS, estimator="BMA")),
               predict(hald.ZS, estimator="BMA")$fit)
  expect_equal(as.vector(fitted(hald.ZS, estimator="HPM")),
               as.vector(predict(hald.ZS, estimator="HPM")$fit))
  expect_equal(as.vector(fitted(hald.ZS, estimator="BPM")),
               as.vector(predict(hald.ZS, estimator="BPM")$fit))
  expect_equal(as.vector(fitted(hald.ZS, estimator="MPM")),
               as.vector(predict(hald.ZS, estimator="MPM")$fit))
})

test_that("GLM logit probne0", {
  data(Pima.tr, package="MASS")
  pima.BAS = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
                     method="BAS",
                     betaprior=bic.prior(), family=binomial(),
                     modelprior=uniform())
  pima.det = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
                     method="deterministic",
                     betaprior=bic.prior(), family=binomial(),
                     modelprior=uniform())
  expect_equal(pima.BAS$probne0, pima.BAS$probne0)

})
