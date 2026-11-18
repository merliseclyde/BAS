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

  # fixed with GitHub issue #52
  expect_no_error(predict(pima_gprior,
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

# related to github issue #93 for default options listed in help
test_that("predict.bas.lm links", {
  data("Hald")
  hald_gprior <- bas.lm(Y ~ ., data = Hald, alpha = 13,
                        prior = "g-prior")
  hald_pred_link <- predict(hald_gprior, newdata = Hald, type = "link",
                           se.fit = TRUE)
  
  hald_pred_response <- predict(hald_gprior, newdata = Hald, type = "response",
                                se.fit = TRUE)
  expect_equal(hald_pred_link$fit, hald_pred_response$fit)
})
  
# GitHub issue #70: exact example from the issue
test_that("predict with se.fit = TRUE and BMA when n < p", {
  data("bodyfat")
  for (method in c("BAS", "deterministic")) {
    bas_mod <- bas.lm(Bodyfat ~ ., data = bodyfat[1:14, ], method = method)
    # models with shrinkage > 1 under the ZS-null Laplace approximation
    # produce NaN for model specific SEs, but have negligible posterior
    # probability
    expect_no_error(pred <- suppressWarnings(
      predict(bas_mod,
        newdata = bodyfat[15:20, ], se.fit = T,
        estimator = "BMA"
      )
    ))
    expect_length(pred$fit, 6)
    expect_true(all(is.finite(pred$se.bma.fit)))
    expect_true(all(is.finite(pred$se.bma.pred)))
    expect_true(all(pred$se.bma.pred > pred$se.bma.fit))
    ci <- confint(pred)
    expect_is(ci, "confint.bas")
    expect_true(all(ci[, 1] < ci[, 2]))
  }
})

# GitHub issue #95: exact example from the issue
test_that("predict.basglm with se.fit = TRUE", {
  data(Pima.tr, package = "MASS")

  # enumeration  with default method="BAS"
  pima.cch <- bas.glm(type ~ .,
    data = Pima.tr, n.models = 2^7,
    method = "BAS",
    betaprior = CCH(a = 1, b = 532 / 2, s = 0), family = binomial(),
    modelprior = beta.binomial(1, 1)
  )

  expect_no_error(pred <- predict(pima.cch,
    newdata = Pima.tr[1:10, ],
    type = "response", se.fit = TRUE
  ))
  expect_length(pred$fit, 10)
  expect_length(pred$se.bma.fit, 10)
  expect_equal(dim(pred$se.fit), c(pima.cch$n.models, 10))
  expect_true(all(is.finite(pred$se.fit)) && all(pred$se.fit > 0))
  expect_true(all(is.finite(pred$se.bma.fit)) && all(pred$se.bma.fit > 0))

  # credible intervals for probabilities are on the response scale
  ci <- confint(pred, parm = "mean")
  expect_true(all(ci[, 1:2] >= 0 & ci[, 1:2] <= 1))
  expect_true(all(ci[, 1] < pred$fit & pred$fit < ci[, 2]))
  expect_equivalent(ci[, 3], pred$fit)

  # single models: intervals agree with transforming the intervals for
  # the linear predictor (HPD intervals under BMA are not invariant)
  for (estimator in c("HPM", "BPM")) {
    pred <- predict(pima.cch,
      newdata = Pima.tr[1:10, ], estimator = estimator,
      type = "response", se.fit = TRUE
    )
    pred_link <- predict(pima.cch,
      newdata = Pima.tr[1:10, ], estimator = estimator,
      type = "link", se.fit = TRUE
    )
    ci <- confint(pred, parm = "mean")
    ci_link <- confint(pred_link, parm = "mean")
    expect_true(all(ci[, 1:2] >= 0 & ci[, 1:2] <= 1))
    expect_true(all(ci[, 1] < pred$fit & pred$fit < ci[, 2]))
    expect_equal(unclass(ci)[, 1:2], plogis(unclass(ci_link)[, 1:2]),
                 check.attributes = FALSE)
  }
})
