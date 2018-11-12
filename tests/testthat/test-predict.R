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

test_that("MPM and Heredity", {
   library(MASS)
   data(Pima.tr)
   library(dplyr)
   Pima.tr = mutate(Pima.tr,
                    Z = Pima.tr$age* Pima.tr$bp*5 + rnorm(nrow(Pima.tr)))
   pima_BAS <-  bas.lm(Z ~ (age + bp + npreg)^2,
                        data = Pima.tr,
                        method = "BAS",
                        prior = "BIC",
                        update = NULL,
                        modelprior = uniform(),
                        force.heredity = TRUE, pivot=TRUE)
    pima_BAS_nc <-  bas.lm(Z ~ (age + bp + npreg)^2,
                           data = Pima.tr, method = "deterministic",
                           prior = "BIC",
                           update = NULL,
                           modelprior = uniform(),
                           force.heredity = FALSE, pivot=TRUE)
    pima_BAS_no <- force.heredity.bas(pima_BAS_nc)
    skip('env issues')  # FIXME
    predict(pima_BAS_nc, newdata=Pima.tr, estimator="MPM", se.fit=TRUE)
    predict(pima_BAS, estimator="MPM", se.fit=TRUE)
})

# question about env
test_that("LM and env", {
  data(Pima.tr, package="MASS")
  library(dplyr)
  Pima.tr = mutate(Pima.tr,
                   Z = Pima.tr$age* Pima.tr$bp*5 + rnorm(nrow(Pima.tr)))
  pima_lm <-  lm(Z ~ (age + bp + npreg)^2,
                      data = Pima.tr)
  skip('issue with env')
  expect_error(lm(eval(pima_lm$call$formula, parent.frame()),
     data=eval(pima_lm$call$data, parent.frame(), parent.frame())))
#  lm( eval(pima_lm$call, env=parent.frame()))
#  call.org = pima_lm$call
#  call.org$formula = as.numeric(type) ~ glu + bp
#  tmp = lm(eval(call.org, parent.frame()))
#  tmp2 = lm(eval(tmp$call))

})

# FIXME
test_that("predict.bas.glm", {

  data(Pima.tr, package="MASS")
  data(Pima.te, package="MASS")
  pima_gprior <- bas.glm(type ~ . data=Pima.tr,
                         betaprior = g.prior(g = as.numeric(nrow(Pima.tr))),
                         family=binomial())


  # should not error
  expect_error(predict(pima_gprior,
                       estimator = "HPM",
                       se.fit = TRUE))
  expect_error(plot(confint(pima_pred, parm = "mean")))
#  should not error
  expect_error( predict(hald_gprior, newdata=Pima.te, estimator = "HPM",
                       se.fit = TRUE))
  expect_null(plot(confint(pima_pred)))
})
