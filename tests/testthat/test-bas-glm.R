context("bas.glm")

test_that("bas.glm initprobs" , {
 data(Pima.tr, package="MASS")
 expect_error(bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      initprobs=rep(.4, nrow(Pima.tr)-1),
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  )
  set.seed(1)
  pima_bas1 <- bas.glm(type ~ .,
                     data = Pima.tr, method = "BAS",
                     initprobs=rep(.4, ncol(Pima.tr)-1),
                     betaprior = bic.prior(), family = binomial(),
                     modelprior = uniform())
  set.seed(1)
  pima_bas2 <- bas.glm(type ~ .,
                     data = Pima.tr, method = "BAS",
                     initprobs=c(1,rep(.4, ncol(Pima.tr)-1)),
                     betaprior = bic.prior(), family = binomial(),
                     modelprior = uniform())
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
  set.seed(1)
  pima_bas2 <- bas.glm(type ~ .,
                       data = Pima.tr, method = "BAS",
                       initprobs=c(1.5,rep(.4, ncol(Pima.tr)-1)),
                       betaprior = bic.prior(), family = binomial(),
                       modelprior = uniform())
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
})

test_that("GLM logit", {
  data(Pima.tr, package = "MASS")
  set.seed(1)
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = bic.prior(), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$probne0[-1] > 1))
  set.seed(1)
  pima_BAS2 <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(), family = "binomial",
                      modelprior = uniform()
  )
  expect_equal(pima_BAS$postprobs, pima_BAS2$postprobs)
  expect_error(bas.glm(type ~ .,
                       data = Pima.tr, method = "BAS",
                       betaprior = bic.prior(),
                       family = "gaussian",
                       modelprior = uniform())
  )
  expect_error(bas.glm(type ~ .,
                       data = Pima.tr, method = "BAS",
                       betaprior = bic.prior(),
                       family = "homeless",
                       modelprior = uniform())
  )
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
    method = "deterministic", betaprior = bic.prior(),
    family = binomial(), modelprior = uniform()
  )
  expect_equal(pima_BAS$probne0, pima_det$probne0)
  expect_equal(
    predict(pima_BAS, type = "link")$fit,
    as.vector(fitted(pima_det))
  )
  expect_equal(
    predict(pima_BAS, data = Pima.tr, type = "link", se.fit = TRUE)$se.bma.fit,
    predict(pima_BAS, type = "link", se.fit = TRUE)$se.bma.fit
  )
  expect_equal(
    predict(pima_BAS, data = Pima.tr, type = "response", se.fit = TRUE)$se.bma.fit,
    predict(pima_BAS, type = "response", se.fit = TRUE)$se.bma.fit
  )
  expect_equal(
    predict(pima_BAS, type = "response")$fit,
    fitted(pima_det, type = "response")
  )
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = Jeffreys(), family = binomial(),
    modelprior = tr.beta.binomial(1, 1, 4))
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
    method = "deterministic", betaprior = Jeffreys(),
    family = binomial(), modelprior = tr.beta.binomial(1, 1, 4))
  expect_equal(pima_BAS$probne0, pima_det$probne0)

  pima_BAS <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = Jeffreys(), family = binomial(),
                      modelprior = tr.poisson(2, 4))
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
                      method = "deterministic", betaprior = Jeffreys(),
                      family = binomial(),
                      modelprior = tr.poisson(2, 4))
  expect_equal(pima_BAS$probne0, pima_det$probne0)
})

# FIXME issue #28
test_that("diagnostic plot for glm MCMC", {
  data(Pima.tr, package="MASS")
  pima_MCMC <- bas.glm(type ~ .,
                       data = Pima.tr, MCMC.iterations = 1024,
                       method = "MCMC", betaprior = aic.prior(),
                       family = binomial(),
                       modelprior = tr.poisson(2,5))
  expect_null(diagnostics(pima_MCMC, type = "model"))
  expect_error(diagnostics(pima_MCMC, type = "pip"))
})

# FIXME issue #35
test_that("MCMC+BAS: missing MCMC.iterations and n.models arg", {
  data(Pima.tr, package = "MASS")
  set.seed(1)
  pima_BAS <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(),
                      family = binomial(),
                      modelprior = uniform())
  set.seed(1)
  pima_1 <- bas.glm(type ~ ., data=Pima.tr,
                        method = "MCMC+BAS",
                        betaprior = bic.prior(),
                        family = binomial(),
                        modelprior = uniform())
  set.seed(1)
  pima_2 <- bas.glm(type ~ ., data = Pima.tr,
                        method = "MCMC+BAS",
                        betaprior = bic.prior(),
                        n.models=2^7,
                        MCMC.iterations=10000, #default
                        family = binomial(),
                        modelprior = uniform())
expect_equal(pima_1$probne0, pima_2$probne0)
expect_equal(pima_BAS$probne0, pima_2$probne0)
expect_equal(pima_BAS$n.models, pima_1$n.models)
})

test_that("missing data arg", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  attach(Pima.tr)
  pima_no_data <- bas.glm(type ~ npreg + glu + bp + skin + bmi + ped + age,
                          method = "BAS",
                          betaprior = bic.prior(),
                          family = binomial(),
                          modelprior = uniform())
  expect_equal(pima_BAS$probne0, pima_no_data$probne0)
  })

test_that("poisson regression", {
  data(crabs, package = "glmbb")
  crabs.bas <- bas.glm(satell ~ color * spine * width + weight,
    data = crabs,
    family = poisson(),
    betaprior = EB.local(), modelprior = uniform(),
    method = "MCMC", n.models = 1024, MCMC.iterations = 10000,
    prob.rw = .95
  )
  expect_null(plot(crabs.bas))
  expect_equal(0, sum(crabs.bas$shrinkage > 1))
})

test_that("glm_fit", {
  data(Pima.tr, package = "MASS")
  Y <- as.numeric(Pima.tr$type) - 1
  X <- cbind(1, as.matrix(Pima.tr[, 1:7]))
  pima_new <- bayesglm.fit(X, Y,
    family = binomial(),
    coefprior = bic.prior(n = length(Y))
  )
  pima_orig <- glm(type ~ ., family = binomial(), data = Pima.tr)
  expect_equivalent(pima_new$coefficients, coef(pima_orig))
  expect_equivalent(pima_new$se, summary(pima_orig)$coef[, 2])
  pima_nowts <- bayesglm.fit(X, Y,
    weights = NULL, offset = NULL,
    family = binomial(),
    coefprior = bic.prior(n = length(Y))
  )
  expect_equal(
    pima_new$coefficients,
    pima_nowts$coefficients
  )
})

test_that("robust prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = robust(), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(nrow(Pima.tr), pima_BAS$betaprior$hyper.parameters$n)
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
})

test_that("intrinsic prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = intrinsic(n = nrow(Pima.tr)),
    family = binomial(),
    modelprior = uniform()
  )
  pima_BAS_no_n <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = intrinsic(),
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
  expect_equal(0, sum(pima_BAS_no_n$shrinkage > 1))
  expect_equal(pima_BAS$probne0, pima_BAS_no_n$probne0)
})

test_that("TestBF prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = testBF.prior(g = nrow(Pima.tr)),
    family = binomial(),
    modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
})

test_that("hyper.g.n prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                      betaprior = hyper.g.n(),
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
})

test_that("hyper.g prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                      betaprior = hyper.g(),
                      family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
  expect_error(bas.glm(type ~ ., data = Pima.tr, method = "BAS",
                       betaprior = hyper.g(alpha=2.0),
                       family = binomial(),
                       modelprior = uniform())
  )
})
test_that("g prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = g.prior(g = as.numeric(nrow(Pima.tr))),
    family = binomial(),
    modelprior = uniform()
  )
  expect_equal(0, sum(pima_BAS$shrinkage > 1))
})

test_that("beta.prime prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = beta.prime(n = 200), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(nrow(Pima.tr), pima_BAS$betaprior$hyper.parameters$n)
  expect_error(pima_BAS <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = beta.prime(), family = binomial(),
    modelprior = uniform()
  ))
  #  expect_equal(nrow(Pima.tr), pima_BAS$betaprior$hyper.parameters$n)
})

test_that("cch prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_cch <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = CCH(2, 2), family = binomial(),
    modelprior = uniform()
  )
  pima_TG <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = TG(), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(pima_cch$probne0, pima_TG$probne0)
})

test_that("Tcch prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima_Tcch <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = tCCH(alpha = 2), family = binomial(),
    modelprior = uniform()
  )
  pima_cch <- bas.glm(type ~ ., data = Pima.tr, method = "BAS",
    betaprior = CCH(2, 2), family = binomial(),
    modelprior = uniform()
  )
  expect_equal(pima_cch$probne0, pima_Tcch$probne0)
  pima_tcch <- bas.glm(type ~ ., data = Pima.tr,
                       method = "BAS",
                       betaprior = tCCH(alpha = 1,
                                        beta = 1,
                                        s = .5),
                       family = binomial(),
                       modelprior = uniform()
  )
  expect_equal(0, sum(pima_tcch$shrinkage > 1))
})

test_that("IC.prior", {
  data(Pima.tr, package = "MASS")
  pima_bic <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  pima_ic <- bas.glm(type ~ .,
                      data = Pima.tr, method = "BAS",
                      betaprior = IC.prior(log(nrow(Pima.tr))), family = binomial(),
                      modelprior = uniform())
  expect_equal(pima_bic$probne0, pima_ic$probne0)
  })

test_that("cv.summary", {
  data(Pima.tr, package = "MASS")
  data(Pima.te, package = "MASS")
  pima_bic <- bas.glm(type ~ .,
                      data = Pima.tr, method = "MCMC",
                      MCMC.iterations = 10000,
                      betaprior = bic.prior(), family = binomial(),
                      modelprior = uniform())
  pima_pred <- predict(pima_bic, newdata=Pima.te, type="response")
  expect_equal(TRUE, cv.summary.bas(pima_pred$fit,
                               as.numeric(Pima.te$type)) > 0)
  expect_equal(TRUE, cv.summary.bas(pima_pred$fit,
                                 as.numeric(Pima.te$type),
                                 score="miss-class")
               <1)
  expect_error(cv.summary.bas(pima_pred$fit,
                                    as.numeric(Pima.te$type),
                                    score="percent-explained"))
  expect_error(cv.summary.bas(pima_pred$fit,
                              as.numeric(Pima.te$type[-1]),
                              score="percent-explained"))
})

test_that("include always", {
  data("Pima.tr", package="MASS")
  expect_error(bas.glm(type ~ .,
                        data = Pima.tr, method = "MCMC",
                        include.always = ~ bp,
                        betaprior = g.prior(g=100), family = binomial(),
                        modelprior = beta.binomial(1, 1))
  )
  # issue #34
  # error Error in bas.glm(type ~ ., data = Pima.tr, method = "MCMC", include.always = ~bp,  :
  # object 'prob' not found
})

test_that("Jeffreys & MCMC", {
 data(Pima.tr, package="MASS")
 pima_BAS <-  bas.glm(type ~ .,
                      data = Pima.tr, method = "MCMC",
                      betaprior = Jeffreys(),
                      family = binomial(),
                      modelprior = tr.beta.binomial(1, 1, 4))
# expect_equal(0, sum(pima_BAS$probne0 > 1))
#  issue #33
})
