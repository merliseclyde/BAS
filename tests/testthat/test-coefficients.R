context("coefficients.bas")

test_that("coefficients", {
  data("Hald")
  hald_gprior <- bas.lm(Y ~ ., data = Hald, prior = "ZS-null",
                        modelprior = beta.binomial(1, 1))
  expect_is(coefficients(hald_gprior), "coef.bas")
  expect_length(coefficients(hald_gprior, estimator='MPM'),11)
  expect_equal("HPM", coefficients(hald_gprior, estimator='HPM')$estimator)
  expect_error(coefficients(hald_gprior, estimator='BPM'))
  expect_null(print(coefficients(hald_gprior, estimator="BMA")))
  expect_equal(coefficients(hald_gprior, estimator='HPM')$postmean,
               coefficients(hald_gprior, estimator='BMA', n.models=1)$postmean)
})
