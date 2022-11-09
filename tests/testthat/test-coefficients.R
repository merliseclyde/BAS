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

# GitHub issue #56 and #39

test_that("formula and env issues with MPM",{
  data(UScrime, package = "MASS")
  UScrime <- UScrime[, 1:5]
  
  crime.bic1 <- bas.lm(formula = M ~ So + Ed + Po1 + Po2,
                       data = UScrime,
                       prior = "JZS",
                       initprobs = c(1, 0.5, 0.5, 0.5, 0.5),
                       renormalize = TRUE)
  coef.mpm1 <- coef(crime.bic1, estimator = "MPM")
  
  form <- M ~ So + Ed + Po1 + Po2
  crime.bic2 <- bas.lm(formula = form,
                       data = UScrime,
                       prior = "JZS",
                       initprobs = c(1, 0.5, 0.5, 0.5, 0.5),
                       renormalize = TRUE)
  
  coef.mpm2 = coef(crime.bic2, estimator = "MPM")
  expect_equal(coef.mpm1, coef.mpm2)
}  
)


test_that("env and lm", {
  data(UScrime, package = "MASS")
  UScrime <- UScrime[, 1:5]
  
  mylm = function(object) {
    modelform = as.formula(eval(object$call$formula, parent.frame()))
    environment(modelform) = environment()
    data = eval(object$call$data)
    weights = eval(object$call$weights)
    
    object = lm(formula = modelform,
                data = data,
                weights = weights)
   return(object) }
  
  
  crime.lm1 <- lm(formula = M ~ So + Ed + Po1 + Po2, data = UScrime)  
  tmp1 = mylm(crime.lm1)
  
  form = M ~ So + Ed + Po1 + Po2
  crime.lm2 <- lm(formula = form, data = UScrime)  
  
  tmp = mylm(crime.lm2)
  
  expect_equal(coef(tmp), coef(tmp1))
  
})

