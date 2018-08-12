context("predict.R")

test_that("predict.bas.lm", {
  data("Hald")
  hald_gprior <- bas.lm(Y ~ ., data = Hald, alpha = 13, prior = "g-prior")
  hald_pred <- predict(hald_gprior, newdata = Hald, estimator = "BPM",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred, parm = "mean")))
  hald_pred <- predict(hald_gprior, estimator = "BMA",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred)))
})

test_that("predict.bas.glm", {
  data("Pima.tr", package="MASS")
  data("Pima.te", package="MASS")
  pima_gprior <- bas.glm(type ~ ., data = Pima.tr,
                         betaprior = g.prior(g=as.numeric(nrow(Pima.tr))),
                         family=binomial())
  pima_pred <- predict(pima_gprior,
                       estimator = "HPM",
                       se.fit = FALSE)

  # should not error
  expect_error(predict(pima_gprior,
                       estimator = "HPM",
                       se.fit = TRUE))
#  expect_null(plot(confint(pima_pred, parm = "mean")))
#  should not error
  expect_error( predict(hald_gprior, newdata=Pima.te, estimator = "HPM",
                       se.fit = TRUE))
  expect_error( predict(hald_gprior, data=Pima.te, estimator = "HPM",
                        se.fit = TRUE))
  expect_null(plot(confint(pima_pred)))
})
