context("bas.lm")

test_that("shrinkage is less than or equal to 1", {
  data(Hald)
  hald.bas <- bas.lm(Y ~ .,
    prior = "ZS-null", modelprior = uniform(),
    data = Hald
  )
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas <- bas.lm(Y ~ .,
    prior = "EB-local", modelprior = uniform(),
    data = Hald
  )
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas <- bas.lm(Y ~ .,
    prior = "EB-global", modelprior = uniform(),
    data = Hald
  )
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas <- bas.lm(Y ~ .,
    prior = "hyper-g", modelprior = uniform(),
    data = Hald
  )
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas <- bas.lm(Y ~ .,
    prior = "hyper-g-n", modelprior = uniform(),
    data = Hald
  )
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas <- bas.lm(Y ~ .,
    prior = "hyper-g-laplace", modelprior = uniform(),
    data = Hald
  )
  expect_equal(0, sum(hald.bas$shrinkage > 1))
  hald.bas <- bas.lm(Y ~ .,
    prior = "g-prior", modelprior = uniform(),
    data = Hald
  )
  expect_equal(0, sum(hald.bas$shrinkage > 1))
})

test_that("shrinkage is less than or equal to 1", {
  data(Hald)
  #  hald.bas = bas.lm(Y ~ ., prior="JZS", modelprior=uniform(), data=Hald)
  #  expect_equal(0, sum(hald.bas$shrinkage > 1))
})

test_that("A/BIC: shrinkage is equal to 1", {
  data(Hald)
  hald.BIC <- bas.lm(Y ~ .,
    prior = "BIC", modelprior = uniform(),
    data = Hald
  )
  expect_equal(hald.BIC$n.model, sum(hald.BIC$shrinkage == 1))
  hald.AIC <- bas.lm(Y ~ .,
    prior = "AIC", modelprior = uniform(),
    data = Hald
  )
  expect_equal(hald.AIC$n.model, sum(hald.AIC$shrinkage == 1))
})

test_that("no method", {
  data(Hald)
  expect_error(bas.lm(Y ~ .,
    prior = "garbage",
    modelprior = uniform(), data = Hald
  ))
})

test_that("deterministic, BAS and MCMC+BAS", {
  data(Hald)
  hald.bas <- bas.lm(Y ~ .,
    prior = "BIC",
    modelprior = uniform(), data = Hald
  )
  hald.MCMCbas <- bas.lm(Y ~ .,
    prior = "BIC", method = "MCMC+BAS",
    modelprior = uniform(), data = Hald, MCMC.iterations = 1000
  )
  hald.deterministic <- bas.lm(Y ~ .,
    prior = "BIC",
    method = "deterministic",
    modelprior = uniform(), data = Hald
  )
  expect_equal(hald.bas$probne0, hald.deterministic$probne0)
  expect_equal(hald.bas$probne0, hald.MCMCbas$probne0)
})

test_that("pivot", {
  data(Hald)
  hald.bas <- bas.lm(Y ~ .,
    prior = "BIC",
    modelprior = uniform(), data = Hald
  )
  hald.deterministic <- bas.lm(Y ~ .,
    prior = "BIC",
    method = "deterministic",
    modelprior = uniform(), data = Hald, pivot = TRUE
  )
  expect_equal(hald.bas$probne0, hald.deterministic$probne0)
  expect_length(summary(hald.bas), 60)
  expect_null(print(hald.bas))
})


test_that("pivoting with non-full rank design", {
  set.seed(42)
  dat <- data.frame(Y = rnorm(5), X1 = 1:5, X2 = 1:5, X3 = rnorm(5))

  tmp.bas <- bas.lm(Y ~ .,
    data = dat, prior = "BIC", modelprior = uniform(),
    method = "BAS", pivot = T
  )

  tmp.mcmc <- bas.lm(Y ~ .,
    data = dat, prior = "BIC", modelprior = uniform(),
    method = "MCMC", pivot = T, MCMC.iterations = 10000
  )
  expect_equal(sort(tmp.bas$R2), sort(tmp.mcmc$R2))
})

test_that("prediction versus fitted", {
  data(Hald)
  hald.ZS <- bas.lm(Y ~ .,
    prior = "ZS-null", modelprior = uniform(),
    data = Hald
  )
  expect_equal(
    as.vector(fitted(hald.ZS, estimator = "BMA")),
    predict(hald.ZS, estimator = "BMA", se.fit = TRUE)$fit
  )
  expect_equal(
    as.vector(fitted(hald.ZS, estimator = "HPM")),
    as.vector(predict(hald.ZS, estimator = "HPM", se.fit = TRUE)$fit)
  )
  expect_equal(
    as.vector(fitted(hald.ZS, estimator = "BPM")),
    as.vector(predict(hald.ZS, estimator = "BPM", se.fit = TRUE)$fit)
  )
  expect_equal(
    as.vector(fitted(hald.ZS, estimator = "MPM")),
    as.vector(predict(hald.ZS, estimator = "MPM", se.fit = TRUE)$fit)
  )
})

test_that("methods", {
  data(Hald)
  expect_error(bas.lm(Y ~ .,
    prior = "ZS-null", modelprior = uniform(),
    data = Hald, method = "AMCMC"
  ))
  expect_error(bas.lm(Y ~ .,
    prior = "hyperg/n", modelprior = uniform(),
    data = Hald
  ))
})

test_that("force.heredity", {
  loc <- system.file("testdata", package = "BAS")
  d <- read.csv(paste(loc, "JASP-testdata.csv", sep = "/"))

  simpleFormula <- as.formula("contNormal ~ contGamma + contcor1 +
                               contGamma*contcor1")

  set.seed(1)
  basObj <- bas.lm(simpleFormula,
    data = d,
    alpha = 0.125316,
    prior = "JZS",
    include.always = as.formula("contNormal ~ contcor1"),
    modelprior = beta.binomial(1, 1),
    weights = d$facFifty
  )
  set.seed(1)
  basObj.old <- bas.lm(simpleFormula,
    data = d,
    alpha = 0.125316,
    prior = "JZS",
    include.always = as.formula("contNormal ~ contcor1"),
    modelprior = beta.binomial(),
    weights = d$facFifty, force.heredity = FALSE
  )
  basObj.old <- force.heredity.bas(basObj.old)

  expect_equal(basObj$probne0, basObj.old$probne0)
})


test_that("check non-full rank", {
  loc <- system.file("testdata", package = "BAS")
  d <- read.csv(paste(loc, "JASP-testdata.csv", sep = "/"))

  fullModelFormula <- as.formula("contNormal ~ contGamma*contExpon +
                                  contGamma*contcor1 + contExpon * contcor1")

  expect_warning(bas.lm(fullModelFormula,
    data = d,
    alpha = 0.125316,
    prior = "JZS",
    weights = facFifty, force.heredity = FALSE, pivot = F
  ))
  expect_error(eplogprob(lm(fullModelFormula, data = d)))
  basObj.eplogp <- bas.lm(fullModelFormula,
    data = d,
    alpha = 0.125316, initprobs = "marg-eplogp",
    prior = "JZS", method = "deterministic", pivot = T,
    modelprior = uniform(),
    weights = facFifty, force.heredity = FALSE
  )
  basObj.det <- bas.lm(fullModelFormula,
    data = d,
    alpha = 0.125316,
    modelprior = uniform(),
    prior = "JZS", method = "deterministic", pivot = T,
    weights = facFifty, force.heredity = FALSE
  )
  basObj <- bas.lm(fullModelFormula,
    data = d,
    alpha = 0.125316, modelprior = uniform(),
    prior = "JZS", method = "BAS", pivot = T,
    weights = facFifty, force.heredity = FALSE
  )
  expect_equal(0, sum(is.na(basObj.det$postprobs)))
  expect_equal(basObj.eplogp$probne0, basObj.det$probne0)
  expect_equal(basObj.det$probne0, basObj$probne0)
  expect_equal(basObj.eplogp$probne0, basObj$probne0)

  basObj.EBL <- bas.lm(fullModelFormula,
    data = d,
    alpha = 0.125316, initprobs = "marg-eplogp",
    prior = "EB-local", method = "deterministic", pivot = T,
    modelprior = uniform(),
    weights = facFifty, force.heredity = FALSE
  )
  basObj.up <- update(basObj.eplogp, newprior = "EB-local")
  expect_equal(basObj.EBL$postprobs, basObj.up$postprobs, tolerance = 1e-03)
})

test_that("as_matrix tools", {
  data(Hald)
  Hald.bic <- bas.lm(Y ~ .,
    data = Hald, prior = "BIC",
    initprobs = "eplogp"
  )
  m1 <- which.matrix(Hald.bic$which, Hald.bic$n.vars)
  colnames(m1) <- Hald.bic$namesx
  m2 <- list2matrix.which(Hald.bic)
  expect_equal(m1, m2)
  m3 <- list2matrix.bas(Hald.bic, "which") > 0
  m3[, 1] <- 1
  probne0 <- t(m3) %*% Hald.bic$postprobs
  expect_equal(as.vector(probne0), Hald.bic$probne0)
})

