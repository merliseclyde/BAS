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
                       modelprior = uniform())
  
  expect_equal(pima_bas1$n.models, pima_bas2$n.models)
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
  expect_equal(pima_bas1$logmarg, pima_bas2$logmarg)
  expect_equal(pima_bas1$freq, pima_bas2$freq)
  expect_equal(pima_bas1$probne0, pima_bas2$probne0)

})

test_that("Test MCMC with Growable Vectors when needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  
  set.seed(42)
  crime.mcmc <-  bas.lm(y ~ ., data=UScrime, n.models=2^16, prior="BIC",
                        method="MCMC", MCMC.it = 100000, burnin = 1000)
  
  set.seed(42)
  crime.grow = bas.lm(y ~ ., data=UScrime, prior="BIC", 
                      n.models = crime.mcmc$n.models - 100,
                      method="MCMC_GROWABLE",  
                      MCMC.it = 100000, burnin = 1000)
  
  expect_equal(crime.grow$logmarg, crime.mcmc$logmarg)
  expect_equal(crime.grow$freq, crime.mcmc$freq)
  expect_equal(crime.grow$size, crime.mcmc$size)
  expect_equal(crime.grow$probne0, crime.mcmc$probne0)
})
