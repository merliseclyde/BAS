context("predict.R")

test_that("predict.bas", {
  data("Hald")
  hald_gprior <- bas.lm(Y ~ ., data = Hald, alpha = 13, prior = "g-prior")
  hald_pred <- predict(hald_gprior, newdata = Hald, estimator = "BPM",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred, parm = "mean")))
  hald_pred <- predict(hald_gprior, estimator = "BMA",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred)))
})
