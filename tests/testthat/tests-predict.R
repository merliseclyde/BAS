context("predict.R")

test_that("predict.bas", {
  data("Hald")
  hald_gprior <- bas.lm(Y ~ ., data = Hald, alpha = 13, prior = "g-prior")
  hald_pred <- predict(hald_gprior, newdata = Hald, estimator = "BPM",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred, parm = "mean")))
  expect_null(plot(confint(hald_pred)))
  hald_pred <- predict(hald_gprior, estimator = "BMA", predict = TRUE,
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred)))
})
