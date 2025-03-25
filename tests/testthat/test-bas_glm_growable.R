context("bas.glm Growable vectors")

test_that("Test MCMC with Growable Vectors when not needed", {
  data(Pima.tr, package="MASS")
  # issue #91 implement growable vectors in MCMC_GROWABLE
  # 
  set.seed(1)
  pima_bas2 <- bas.glm(type ~ .,
                       data = Pima.tr, 
                       method="MCMC",  MCMC.it = 10000, burnin = 1000,
                       initprobs=c(1,rep(.4, ncol(Pima.tr)-1)),
                       betaprior = bic.prior(), family = binomial(),
                       modelprior = uniform())
  set.seed(1)
  pima_bas1 <- bas.glm(type ~ .,
                       data = Pima.tr, method="MCMC_GROWABLE",  
                       MCMC.it = 10000, burnin = 1000,  n.models = pima_bas2$n.models,
                       initprobs=rep(.4, ncol(Pima.tr)-1),
                       betaprior = bic.prior(), family = binomial(),
                       modelprior = uniform(), renormalize = FALSE)
  
  expect_equal(pima_bas1$n.models, pima_bas2$n.models)
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
  expect_equal(pima_bas1$logmarg, pima_bas2$logmarg)
  expect_equal(pima_bas1$freq, pima_bas2$freq)
  expect_equal(pima_bas1$probne0, pima_bas2$probne0)

})

test_that("Test MCMC with Growable Vectors when needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE

})
