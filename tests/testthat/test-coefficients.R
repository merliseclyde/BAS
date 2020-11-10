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

  coef_hald <- coef(hald_gprior)
  confint_hald <- confint(coef_hald)
  expect_is(confint_hald, "confint.bas")
  expect_is(confint(coef_hald, approx=FALSE, nsim=5000), "confint.bas")
  expect_length(confint(coef_hald, parm="X4"), 3)
  expect_null(plot(confint(coef_hald, parm=2:5)))
  expect_null(plot(confint(coef(hald_gprior, estimator='HPM'))))
  expect_null(plot(confint(coef(hald_gprior, estimator='HPM')),
                   horizontal = TRUE))
  expect_error(confint(coef(hald_gprior), nsim=1))
})

test_that("plot posterior coefficients", {
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  crime_bic <- bas.lm(y ~ ., data=UScrime, n.models=2^10, prior="BIC")
  crime_coef <- coef(crime_bic)
  expect_null(plot(crime_coef, subset=2:4, ask=TRUE))
  expect_null(plot(crime_coef, ask=FALSE))
})


test_that("BPM coefficients", {
  # unit test for merliseclyde/BAS#49
  data("Hald")
  hald.bas =  bas.lm(Y~ .,
                     data = Hald,
                     alpha = 1,
                     modelprior = uniform(),
                     prior = "BIC", n.models = 2^4)

  hald.pred <- predict(hald.bas, estimator = "BPM")
  BPM = as.vector(which.matrix(hald.bas$which[hald.pred$best],
                               hald.bas$n.vars))

  hald.BPM = bas.lm(Y ~ ., data = Hald,
                    alpha = 1,
                    prior = "BIC",
                    modelprior = uniform(),
                    bestmodel = BPM, n.models = 1)
  hald.coef <- coef(hald.BPM)

 expect_equal(hald.coef$postmean[BPM == 1], hald.bas$mle[[hald.pred$best]])
})
