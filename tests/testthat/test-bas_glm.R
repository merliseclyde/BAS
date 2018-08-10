context("bas.glm")

test_that("GLM logit", {
  data(Pima.tr, package = "MASS")
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = bic.prior(), family = binomial(),
    modelprior = uniform()
  )
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
    method = "deterministic", betaprior = bic.prior(),
    family = binomial(), modelprior = uniform()
  )
  pima_MCBAS <- bas.glm(type ~ ., data = Pima.tr,
    method = "MCMC+BAS", betaprior = bic.prior(),
    family = binomial(), modelprior = uniform()
  )
  pima_MCMC <- bas.glm(type ~ .,
    data = Pima.tr, MCMC.iterations = 1024,
    method = "MCMC", betaprior = aic.prior(),
    family = binomial(),
    modelprior = tr.poisson(ncol(Pima.tr) - 1, 7)
  )
  expect_null(diagnostics(pima_MCMC, type = "model"))
  expect_error(diagnostics(pima_MCMC, type = "pip"))
  expect_equal(pima_BAS$probne0, pima_det$probne0)
  expect_equal(pima_BAS$probne0, pima_MCBAS$probne0)
  expect_equal(
    predict(pima_BAS, type = "link")$fit,
    as.vector(fitted(pima_det))
  )
  expect_equal(
    predict(pima_BAS, data = Pima.tr, type = "link", se.fit = TRUE)$se.bma.fit,
    predict(pima_BAS, type = "link", se.fit = TRUE)$se.bma.fit
  )
  expect_equal(
    predict(pima_BAS, type = "response")$fit,
    fitted(pima_det, type = "response")
  )
  pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = Jeffreys(), family = binomial(),
    modelprior = tr.beta.binomial(1, 1, 4)
  )
  pima_det <- bas.glm(type ~ ., data = Pima.tr,
    method = "deterministic", betaprior = Jeffreys(),
    family = binomial(), modelprior = tr.beta.binomial(1, 1, 4)
  )
  expect_equal(pima_BAS$probne0, pima_det$probne0)
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
  expect_error(pima_BAS <- bas.glm(type ~ .,
    data = Pima.tr, method = "BAS",
    betaprior = intrinsic(n = nrow(Pima.tr)),
    family = binomial(),
    modelprior = uniform()
  ))
  #  expect_equal(0, sum(pima_BAS$shrinkage > 1))
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
  expect_equal(pima_cch$probne0, pima_Tcch$probne0) <
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
