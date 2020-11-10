context("predict.R")

test_that("predict.bas.lm", {
  data("Hald")
  hald_gprior <- bas.lm(Y ~ ., data = Hald, alpha = 13,
                        prior = "g-prior")
  hald_pred_new <- predict(hald_gprior, newdata = Hald, estimator = "BPM",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred_new, parm = "mean")))
  hald_pred <- predict(hald_gprior, estimator = "BMA",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred)))
  hald_pred_new <- predict(hald_gprior, newdata = Hald,
                           estimator = "BMA",
                           se.fit = TRUE)
  expect_equal(hald_pred_new$fit, hald_pred$fit)
  hald_pred <- predict(hald_gprior, estimator = "HPM",
                       se.fit = TRUE)
  expect_null(plot(confint(hald_pred)))
  hald_pred_new <- predict(hald_gprior, newdata = Hald,
                           estimator = "HPM",
                           se.fit = TRUE)
  expect_equal(hald_pred_new$fit, hald_pred$fit)
  hald_pred_new <- predict(hald_gprior, newdata = Hald[1,],
                           estimator = "HPM",
                           se.fit = TRUE)
  expect_equivalent(hald_pred_new$fit, hald_pred$fit[1])
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
  pima_top <-  predict(pima_gprior,
                       estimator = "BMA", top=1,
                       se.fit = TRUE)

 expect_equal(pima_pred$fit, pima_top$fit, check.attributes = FALSE)

})



#Fixed Issue #51
test_that("MPM and predict glm", {
  data("Pima.tr", package="MASS")
  data("Pima.te", package="MASS")
  pima_gprior <- bas.glm(type ~ ., data = Pima.tr,
                         betaprior = g.prior(g=as.numeric(nrow(Pima.tr))),
                         family=binomial())
  pima_MPM = extract_MPM(pima_gprior)

  expect_equal(predict(pima_gprior, estimator = "MPM", se.fit = FALSE)$fit,
               predict(pima_MPM, se.fit = FALSE)$fit,
               check.attributes = FALSE)


  pima_pred <- predict(pima_gprior,
                       estimator = "MPM", type = "link",
                       se.fit = FALSE)
  pima_fit <-  fitted(pima_gprior,
                       estimator = "MPM")

  expect_equal(pima_pred$fit, pima_fit, check.attributes = FALSE)

})


# Issue #52 SE's are incorrect for glms and weighted regression
test_that("se.fit.glm", {
  data("Pima.tr", package="MASS")
  data("Pima.te", package="MASS")

pima.bic = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
                             method="BAS",
                             betaprior=bic.prior(n=200), family=binomial(),
                             modelprior=beta.binomial(1,1))

fit.bic = predict(pima.bic,  se.fit = TRUE, top=1, type="link", estimator="HPM")
pred.bic = predict(pima.bic, newdata=Pima.te, se.fit = TRUE, top=1, type="link")

form = paste("type ~ ",
             paste0((pred.bic$best.vars[pred.bic$bestmodel[[1]] + 1])[- 1],
                     collapse = "+"))

pima.glm = glm(form, data=Pima.tr, family=binomial())
fit.glm = predict(pima.glm,  se.fit=TRUE, type='link')
pred.glm = predict(pima.glm, newdata=Pima.te, se.fit=TRUE, type='link')

expect_true(all.equal(fit.glm$fit, fit.bic$fit, check.attributes = FALSE))


# issue #50 in github regarding se.fit failing; debugging indicates se.fit is
# incorrect
# Should be expect_equal

expect_equal(fit.glm$se.fit,  fit.bic$se.fit, check.attributes = FALSE)
expect_equal(pred.glm$se.fit, pred.bic$se.fit, check.attributes = FALSE)


})

# Added feature issue #53
test_that("MPM and predict in lm", {
  data(Hald, package="BAS")
  hald_bic =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="BIC",
                     modelprior = uniform())

  hald_MPM = extract_MPM(hald_bic)
  expect_equal(predict(hald_bic, estimator = "MPM")$fit,
               predict(hald_MPM)$fit, check.attributes = FALSE)
  expect_equal(predict(hald_bic, estimator = "MPM")$fit,
               predict(hald_bic, estimator = "MPMold")$fit,
               check.attributes = FALSE)
})
