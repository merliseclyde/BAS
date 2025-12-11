context("bas.glm Growable vectors")

# skip("skip test of bas.glm with Growable Vectors when not needed")
test_that("Test BAS with Growable Vectors when not needed", {
  data(Pima.tr, package="MASS")
  # issue #91 implement growable vectors in MCMC_GROWABLE
  # 
  set.seed(1)
  pima_bas2 <- bas.glm(type ~ .,
                       data = Pima.tr, 
                       method="BAS",  
                       initprobs=c(1,rep(.4, ncol(Pima.tr)-1)),
                       betaprior = bic.prior(), family = binomial(),
                       modelprior = uniform(), GROW = TRUE)
  set.seed(1)
 pima_bas1 <- bas.glm(type ~ .,
                       data = Pima.tr, method="BAS",  
                       initprobs=rep(.4, ncol(Pima.tr)-1),
                       betaprior = bic.prior(), family = binomial(),
                       modelprior = uniform(), GROW = FALSE)
  expect_equal(pima_bas1$n.models, pima_bas2$n.models)
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
  expect_equal(pima_bas1$logmarg, pima_bas2$logmarg)
  expect_equal(pima_bas1$probne0, pima_bas2$probne0)
  
})

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
                       modelprior = uniform(), renormalize = FALSE, GROW = FALSE)
  set.seed(1)
  pima_bas1 <- bas.glm(type ~ .,
                       data = Pima.tr, method="MCMC",  
                       MCMC.it = 10000, burnin = 1000,  n.models = pima_bas2$n.models,
                       initprobs=rep(.4, ncol(Pima.tr)-1),
                       betaprior = bic.prior(), family = binomial(),
                       modelprior = uniform(), renormalize = FALSE, GROW = TRUE)
  
  expect_equal(pima_bas1$n.models, pima_bas2$n.models)
  expect_equal(pima_bas1$postprobs, pima_bas2$postprobs)
  expect_equal(pima_bas1$logmarg, pima_bas2$logmarg)
  expect_equal(pima_bas1$freq, pima_bas2$freq)
  expect_equal(pima_bas1$probne0, pima_bas2$probne0)

})


test_that("Test MCMC with Growable Vectors when needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE

})

test_that("poisson regression GROWABLE", {
  data(crabs, package = "glmbb")
  set.seed(1)
  crabs.MCMCbas <- bas.glm(satell ~ color * spine * width + weight,
                           data = crabs,
                           family = poisson(),
                           betaprior = EB.local(), modelprior = uniform(),
                           method = "MCMC+BAS", 
                           n.models = 1024, burnin = 1000, MCMC.iterations = 10000,
                           prob.rw = .95
  )
  
  set.seed(1)
  crabs.MCMC <- bas.glm(satell ~ color * spine * width + weight,
                        data = crabs,
                        family = poisson(),
                        betaprior = EB.local(), modelprior = uniform(),
                        method = "MCMC", 
                        burnin = 1000, MCMC.iterations = 10000,
                        prob.rw = .95, GROW = FALSE
  )
  set.seed(1)
  crabs.MCMCbas <- bas.glm(satell ~ color * spine * width + weight,
                           data = crabs,
                           family = poisson(),
                           betaprior = EB.local(), modelprior = uniform(),
                           method = "MCMC+BAS", 
                           n.models = crabs.MCMC$n.models, burnin = 1000+10000+1, MCMC.iterations = 0,
                           prob.rw = .95
  )
  set.seed(1)
  crabs.MCMCGrow <- bas.glm(satell ~ color * spine * width + weight,
                            data = crabs,
                            family = poisson(),
                            betaprior = EB.local(), modelprior = uniform(),
                            method = "MCMC", 
                            n.models = 2999, burnin = 1000, MCMC.iterations = 10000,
                            prob.rw = .95, GROW = TRUE
  )
  
  expect_null(plot(crabs.MCMCGrow))
  expect_equal(0, sum(crabs.MCMC$shrinkage > 1))
  expect_equal(0, sum(crabs.MCMCGrow$shrinkage > 1))
  expect_equal(crabs.MCMC$freq, crabs.MCMCGrow$freq)
  expect_equal(crabs.MCMC$probne0, crabs.MCMCGrow$probne0)
  expect_equal(crabs.MCMC$postprobs.MCMC,crabs.MCMCGrow$postprobs)
  expect_equal(crabs.MCMC$probne0.MCMC,crabs.MCMCGrow$probne0.MCMC)
  expect_equal(crabs.MCMCbas$logmarg, crabs.MCMCGrow$logmarg)
  # currently not equal but not used  
  #  expect_equal(crabs.MCMC$freq, crabs.MCMCbas$freq)  
})

