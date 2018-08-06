context("bas.lm")

test_that("shrinkage is less than or equal to 1", {
  data(Hald)
  hald.bas = bas.lm(Y ~ ., prior="JZS", modelprior=uniform(), data=Hald)
  expect_equal(0, sum(hald.bas$shrinkage > 1))
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
