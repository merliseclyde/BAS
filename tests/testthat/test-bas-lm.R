context("bas.lm")

test_that("initprobs out of range", {
  data(Hald)
  bas_hald1 <- bas.lm(Y ~ ., data=Hald, prior="BIC",
                      method="BAS",
                     initprobs = c(-.4, .3, 1.5, .8))
  bas_hald2 <- bas.lm(Y ~ ., data=Hald, prior="BIC",
                      initprobs = c(1, -.4, .3, 1.0, .8))
  bas_hald3 <- bas.lm(Y ~ X1 + X2 + X3 + X4,
                      include.always=~1 + X3,
                      initprobs =  c(1, -.4, .3, 1.0, .8),
                      data=Hald, prior="BIC")
  expect_equal(bas_hald2$probne0, bas_hald2$probne0)
  expect_equal(bas_hald2$probne0, bas_hald3$probne0)
  bas_hald1 <- bas.lm(Y ~ ., data=Hald, prior="BIC",
                      method="MCMC+BAS",
                      initprobs = c(-.4, .3, 1.5, .8))
  bas_hald3 <- bas.lm(Y ~ X1 + X2 + X3 + X4,
                      method="MCMC+BAS",
                      include.always=~1 + X3,
                      initprobs =  c(1, -.4, .3, 1.0, .8),
                      data=Hald, prior="BIC")
  expect_equal(bas_hald2$probne0, bas_hald2$probne0)
  expect_equal(bas_hald2$probne0, bas_hald3$probne0)
  bas_hald1 <- bas.lm(Y ~ ., data=Hald, prior="BIC",
                      method="deterministic",
                      initprobs = c(-.4, .3, 1.5, .8))
  bas_hald3 <- bas.lm(Y ~ X1 + X2 + X3 + X4,
                      method="deterministic",
                      include.always=~1 + X3,
                      initprobs =  c(1, -.4, .3, 1.0, .8),
                      data=Hald, prior="BIC")
  expect_equal(bas_hald2$probne0, bas_hald2$probne0)
  expect_equal(bas_hald2$probne0, bas_hald3$probne0)
})

test_that("shrinkage is less than or equal to 1", {
  data(Hald)
  hald_bas <- bas.lm(Y ~ ., prior = "ZS-null",
                     modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "EB-local", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "EB-global", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "hyper-g", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "hyper-g-n", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "hyper-g-laplace", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "g-prior", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "AIC", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "BIC", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  hald_bas <- bas.lm(Y ~ ., prior = "ZS-full", modelprior = uniform(), data = Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
})

test_that("shrinkage is less than or equal to 1", {
  data(Hald)
  hald_bas = bas.lm(Y ~ ., prior="JZS", modelprior=uniform(), data=Hald)
  expect_equal(0, sum(hald_bas$shrinkage > 1))
  # GitHub  issue # 27
})

test_that("A/BIC: shrinkage is equal to 1", {
  data(Hald)
  hald_BIC <- bas.lm(Y ~ ., prior = "BIC", modelprior = uniform(), data = Hald)
  expect_equal(hald_BIC$n.model, sum(hald_BIC$shrinkage == 1))
  hald_AIC <- bas.lm(Y ~ ., prior = "AIC", modelprior = uniform(), data = Hald)
  expect_equal(hald_AIC$n.model, sum(hald_AIC$shrinkage == 1))
})

test_that("prior not implemented", {
  data(Hald)
  expect_error(bas.lm(Y ~ .,
    prior = "garbage",
    modelprior = uniform(), data = Hald
  ))
})

test_that("inits wrong length", {
  data(Hald)
  expect_error(bas.lm(Y ~ .,
                      prior = "BIC",
                      initprobs = c(1, rep(.4, 8)),
                      modelprior = uniform(), data = Hald
  ))
})

test_that("bernoulli prob wrong length", {
  data(Hald)
  expect_error(bas.lm(Y ~ .,
                      prior = "BIC",
                      modelprior = Bernoulli(probs=rep(.5, 8)), data = Hald
  ))
})

test_that("deterministic, BAS and MCMC+BAS", {
  data(Hald)
  hald_bas <- bas.lm(Y ~ .,
    prior = "BIC",
    modelprior = uniform(), data = Hald
  )

  hald_MCMCbas <- bas.lm(Y ~ .,
                         prior = "BIC", method = "MCMC+BAS",
                         modelprior = uniform(), data = Hald)
                        #MCMC.iterations = 1000)
  hald_deterministic <- bas.lm(Y ~ .,
    prior = "BIC",
    method = "deterministic",
    modelprior = uniform(), data = Hald
  )
  expect_equal(hald_bas$probne0, hald_deterministic$probne0)
  expect_equal(hald_bas$probne0, hald_MCMCbas$probne0)
})

test_that("pivot", {
  data(Hald)
  hald_bas <- bas.lm(Y ~ .,
    prior = "BIC",
    modelprior = uniform(), data = Hald
  )
  hald_deterministic <- bas.lm(Y ~ .,
    prior = "BIC",
    method = "deterministic",
    modelprior = uniform(), data = Hald, pivot = TRUE
  )
  expect_equal(hald_bas$probne0, hald_deterministic$probne0)
})


test_that("pivoting with non-full rank design", {
  set.seed(42)
  dat <- data.frame(Y = rnorm(5), X1 = 1:5, X2 = 1:5, X3 = rnorm(5))

  tmp.bas <- bas.lm(Y ~ ., data = dat, prior = "BIC", modelprior = uniform(), method = "BAS", pivot = T)

  tmp.mcmc <- bas.lm(Y ~ ., data = dat, prior = "BIC", modelprior = uniform(), method = "MCMC", pivot = T, MCMC.iterations = 10000)
  expect_equal(sort(tmp.bas$R2), sort(tmp.mcmc$R2))
})

test_that("prediction versus fitted", {
  data(Hald)
  hald_ZS <- bas.lm(Y ~ .,
    prior = "ZS-null", modelprior = uniform(),
    data = Hald
  )
  expect_equal(
    as.vector(fitted(hald_ZS, estimator = "BMA")),
    predict(hald_ZS, estimator = "BMA", se.fit = TRUE)$fit
  )
  expect_equal(
    as.vector(fitted(hald_ZS, estimator = "HPM")),
    as.vector(predict(hald_ZS, estimator = "HPM", se.fit = TRUE)$fit)
  )
  expect_equal(
    as.vector(fitted(hald_ZS, estimator = "BPM")),
    as.vector(predict(hald_ZS, estimator = "BPM", se.fit = TRUE)$fit)
  )
  expect_equal(
    as.vector(fitted(hald_ZS, estimator = "MPM")),
    as.vector(predict(hald_ZS, estimator = "MPM", se.fit = TRUE)$fit)
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

test_that("big model space", {
  data("protein")
  expect_error(bas.lm(prot.act1 ~ (buf + pH + NaCl + con + ra +
                                     det + MgCl2 + temp)^2,
                      data = protein, n.models = 2^26,
                      method = "BAS", force.heredity = FALSE))
})

test_that("force.heredity", {
  # based on bug #26
  loc <- system.file("testdata", package = "BAS")
  d <- read.csv(paste(loc, "JASP-testdata.csv", sep = "/"))

  simpleFormula <- as.formula("contNormal ~ contGamma + contcor1 + contGamma * contcor1 ")

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

test_that("interactions", {
  data(Hald)
  bas_hald <- bas.lm(Y ~ .^2, data=Hald, method="MCMC",
                     force.heredity = TRUE)
  expect_equal(sum(bas_hald$postprobs >0),
               bas_hald$n.models)
  bas_hald <- bas.lm(Y ~ .^2, data=Hald, method="MCMC",
                     MCMC.iterations = 10000,
                     force.heredity = FALSE)
  bas_hald <- force.heredity.bas(bas_hald)
  expect_equal(sum(bas_hald$postprobs >0),
               bas_hald$n.models)
 })

test_that("check non-full rank", {
  loc <- system.file("testdata", package = "BAS")
  d <- read.csv(paste(loc, "JASP-testdata.csv", sep = "/"))

  fullModelFormula <- as.formula("contNormal ~  contGamma * contExpon + contGamma * contcor1 + contExpon * contcor1")

  bas.lm(fullModelFormula,
         data = d,
         alpha = 0.125316,
         prior = "JZS",
         weights = facFifty, force.heredity = FALSE, pivot = F
  )
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
  expect_equal(basObj.EBL$postprobs, basObj.up$postprobs,
               tolerance=.001)
  # follow up
})

# add more testing for update.bas
#
test_that("as.matrix tools", {
  data(Hald)
  hald_bic <- bas.lm(Y ~ .,
    data = Hald, prior = "BIC",
    initprobs = "eplogp"
  )
  m1 <- which.matrix(hald_bic$which, hald_bic$n.vars)
  colnames(m1) <- hald_bic$namesx
  m2 <- list2matrix.which(hald_bic)
  expect_equal(m1, m2)
  m3 <- list2matrix.bas(hald_bic, "which") > 0
  m3[, 1] <- 1
  probne0 <- t(m3) %*% hald_bic$postprobs
  expect_equal(as.vector(probne0), hald_bic$probne0)
})
