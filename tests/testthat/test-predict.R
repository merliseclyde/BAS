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
  hald_pred_new <- predict(hald_gprior, newdata = Hald, estimator = "BPM",
                           se.fit = FALSE)
  expect_warning(confint(hald_pred_new ))
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
  expect_warning(confint(pima_pred))

  # should not error
  expect_error(predict(pima_gprior,
                       estimator = "HPM",
                       se.fit = TRUE))
#  expect_null(plot(confint(pima_pred, parm = "mean")))
#  should not error
  expect_error( predict(hald_gprior, newdata=Pima.te, estimator = "HPM",
                       se.fit = TRUE))
  #expect_null(plot(confint(pima_pred)))
})

# GitHub issue  #68
# a model with one predictor variable and se.fit=TRUE 
test_that("se.fit with 1 variable", {
  data("Hald")
  hald.gprior =  bas.lm(Y ~ X2, data=Hald, alpha=13, prior="g-prior")
  expect_no_error(predict(hald.gprior, 
                       newdata=Hald, estimator="BPM", se.fit=TRUE))
  
  expect_no_error(predict(hald.gprior, 
                       newdata=Hald, estimator="HPM", se.fit=TRUE))
  
  expect_no_error(predict(hald.gprior, 
                       newdata=Hald, estimator="MPM", se.fit=TRUE))
  
})

# GitHub issue  #70 and #74
test_that("bas.lm using truncated priors includes models with prior prob 0", {
  data("bodyfat")
  bas_mod <- bas.lm(Bodyfat ~.,data = bodyfat[1:14,], method = 'BAS', modelprior = tr.poisson(2, 3))
  
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'HPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'MPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BMA'))
  
  bas_mod <- bas.lm(Bodyfat ~.,data = bodyfat[1:14,], method = 'deterministic', modelprior = tr.poisson(2, 3))
  
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'HPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'MPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BMA'))
  
  bas_mod <- bas.lm(Bodyfat ~.,data = bodyfat[1:14,], method = 'MCMC', modelprior = tr.poisson(2, 3))
  
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'HPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'MPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BMA'))
  
  bas_mod <- bas.lm(Bodyfat ~.,data = bodyfat[1:14,], method = 'MCMC+BAS', modelprior = tr.poisson(2, 3))
  
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'HPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'MPM'))
  expect_no_error(predict(bas_mod,newdata = bodyfat[15:20,], se.fit = T, estimator = 'BMA'))
  


})

# github issue #74
test_that("bas.glm using truncated priors includes models with prior prob 0", {
  data("Pima.tr", package="MASS")
  data("Pima.te", package="MASS")
  bas_mod <- bas.glm(type ~ ., data = Pima.tr, subset = 1:5, method = 'BAS',
                         modelprior = tr.poisson(2,2),
                         betaprior = g.prior(g=as.numeric(nrow(Pima.tr))),
                         family=binomial())
  
  expect_no_error(sum(bas_mod$postprobs == 1.0))
  
# github issue ??
#  expect_no_error(predict(bas_mod,newdata = Pima.tr[15:20,], se.fit = T, estimator = 'HPM'))
#  expect_no_error(predict(bas_mod,newdata = Pima.tr[15:20,], se.fit = T, estimator = 'BPM'))
#  expect_no_error(predict(bas_mod,newdata = Pima.tr[15:20,], se.fit = T, estimator = 'MPM'))
#  expect_no_error(predict(bas_mod,newdata = Pima.tr[15:20,], se.fit = T, estimator = 'BMA'))
  
  #expect_null(plot(confint(pima_pred)))
})