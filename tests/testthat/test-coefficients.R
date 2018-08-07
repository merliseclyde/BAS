context("coefficients.bas")

test_that("coefficients", {
  data("Hald")
  hald.gprior <- bas.lm(Y ~ ., data = Hald, prior = "ZS-null",
                        modelprior = beta.binomial(1,1))
  expect_is(coefficients(hald.gprior), "coef.bas")
})
