context("bas.glm")

test_that("GLM logit", {
  data(Pima.tr, package = "MASS")
  pima.BAS <- bas.glm(type ~ .,
    data = Pima.tr, n.models = 2^7, method = "BAS",
    betaprior = bic.prior(), family = binomial(),
    modelprior = uniform()
  )
  pima.det <- bas.glm(type ~ .,
    data = Pima.tr, n.models = 2^7,
    method = "deterministic", betaprior = bic.prior(),
    family = binomial(), modelprior = uniform()
  )
  pima.MCBAS <- bas.glm(type ~ .,
    data = Pima.tr, n.models = 2^7,
    method = "MCMC+BAS", betaprior = bic.prior(),
    family = binomial(), modelprior = uniform()
  )
  pima.MCMC <- bas.glm(type ~ .,
    data = Pima.tr, MCMC.iterations = 2^10,
    method = "MCMC", betaprior = aic.prior(),
    family = binomial(),
    modelprior = tr.poisson(ncol(Pima.tr) - 1, 7)
  )
  expect_null(diagnostics(pima.MCMC, type = "model"))
  # expect_null(diagnostics(pima.MCMC, type="pip"))
  expect_equal(pima.BAS$probne0, pima.det$probne0)
  expect_equal(pima.BAS$probne0, pima.MCBAS$probne0)
  expect_equal(
    predict(pima.BAS, type = "link")$fit,
    as.vector(fitted(pima.det))
  )
  expect_equal(
    predict(pima.BAS, data = Pima.tr, type = "link", se.fit = TRUE)$se.bma.fit,
    predict(pima.BAS, type = "link", se.fit = TRUE)$se.bma.fit
  )
  expect_equal(
    predict(pima.BAS, type = "response")$fit,
    fitted(pima.det, type = "response")
  )
  pima.BAS <- bas.glm(type ~ .,
    data = Pima.tr, n.models = 2^7, method = "BAS",
    betaprior = Jeffreys(), family = binomial(),
    modelprior = tr.beta.binomial(1, 1, 4)
  )
  pima.det <- bas.glm(type ~ .,
    data = Pima.tr, n.models = 2^7,
    method = "deterministic", betaprior = Jeffreys(),
    family = binomial(), modelprior = tr.beta.binomial(1, 1, 4)
  )
  expect_equal(pima.BAS$probne0, pima.det$probne0)
})

test_that("poisson regression", {
  data(crabs, package = "glmbb")
  crabs.bas <- bas.glm(satell ~ color * spine * width + weight,
    data = crabs,
    family = poisson(),
    betaprior = EB.local(), modelprior = uniform(),
    method = "MCMC", n.models = 2^10, MCMC.iterations = 10000,
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
  pima_nowts <- bayesglm.fit(X, Y, weights=NULL, offset=NULL,
                           family = binomial(),
                           coefprior = bic.prior(n = length(Y)))
  expect_equal(pima_new$coefficients,
                    pima_nowts$coefficients)

})

test_that("robust prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima.BAS <- bas.glm(type ~ .,
                      data = Pima.tr, n.models = 2^7, method = "BAS",
                      betaprior = robust(), family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(nrow(Pima.tr), pima.BAS$betaprior$hyper.parameters$n)
  expect_equal(0, sum(pima.BAS$shrinkage > 1))

  })

test_that("intrinsic prior for GLM", {
  data(Pima.tr, package = "MASS")
  expect_error( pima.BAS <- bas.glm(type ~ .,
                      data = Pima.tr, n.models = 2^7, method = "BAS",
                      betaprior = intrinsic(n = nrow(Pima.tr)), family = binomial(),
                      modelprior = uniform()
  ))
#  expect_equal(0, sum(pima.BAS$shrinkage > 1))
})

test_that("TestBF prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima.BAS <- bas.glm(type ~ .,
                                    data = Pima.tr, n.models = 2^7, method = "BAS",
                                    betaprior = testBF.prior(g = nrow(Pima.tr)), family = binomial(),
                                    modelprior = uniform()
  )
  expect_equal(0, sum(pima.BAS$shrinkage > 1))
})

test_that("beta.prime prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima.BAS <- bas.glm(type ~ .,
                      data = Pima.tr, n.models = 2^7, method = "BAS",
                      betaprior = beta.prime(n=200), family = binomial(),
                      modelprior = uniform()
  )
  expect_equal(nrow(Pima.tr), pima.BAS$betaprior$hyper.parameters$n)
  expect_error(pima.BAS <- bas.glm(type ~ .,
                      data = Pima.tr, n.models = 2^7, method = "BAS",
                      betaprior = beta.prime(), family = binomial(),
                      modelprior = uniform()))
#  expect_equal(nrow(Pima.tr), pima.BAS$betaprior$hyper.parameters$n)
})

test_that("cch prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima.cch <- bas.glm(type ~ .,
                      data = Pima.tr, n.models = 2^7, method = "BAS",
                      betaprior = CCH(2,2), family = binomial(),
                      modelprior = uniform())
  pima.TG <- bas.glm(type ~ .,
                      data = Pima.tr, n.models = 2^7, method = "BAS",
                      betaprior = TG(), family = binomial(),
                      modelprior = uniform())
  expect_equal(pima.cch$probne0, pima.TG$probne0)
})

test_that("Tcch prior for GLM", {
  data(Pima.tr, package = "MASS")
  pima.Tcch <- bas.glm(type ~ .,
                      data = Pima.tr, n.models = 2^7, method = "BAS",
                      betaprior =tCCH(alpha=2), family = binomial(),
                      modelprior = uniform())
  pima.cch <- bas.glm(type ~ .,
                     data = Pima.tr, n.models = 2^7, method = "BAS",
                     betaprior = CCH(2,2), family = binomial(),
                     modelprior = uniform())
  expect_equal(pima.cch$probne0, pima.Tcch$probne0)
})
