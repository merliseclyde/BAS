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

test_that("plot posterior coefficients", {
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  crime_bic <- bas.lm(y ~ ., data=UScrime, n.models=2^10, prior="BIC")
  crime_coef <- coefficients(crime_bic)
  expect_null(plot(crime_coef, subset=2:4, ask=TRUE))
  expect_null(plot(crime_coef, ask=FALSE))
})

