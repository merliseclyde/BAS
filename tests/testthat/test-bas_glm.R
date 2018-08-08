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
})
