context("AMCMC bas.lm")

test_that("AMCMC VS MCMC with no adaptation", {
  data(Hald)
  set.seed(42)
  hald.mcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                     data = Hald, method = "MCMC",  burnin.iteration = 200, MCMC.iterations = 0)
  set.seed(42)
  hald.amcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                      data = Hald, method = "AMCMC", burnin.iteration = 200, MCMC.iterations = 0, 
                      GROW = FALSE)
  expect_equal(hald.amcmc$postprobs.MCMC, hald.mcmc$postprobs.MCMC)
  expect_equal(hald.amcmc$postprobs.MCMC, hald.mcmc$postprobs.MCMC)
  expect_equal(hald.amcmc$n.models, hald.mcmc$n.models)
  expect_equal(hald.amcmc$logmarg, hald.mcmc$logmarg)
  expect_equal(hald.amcmc$freq, hald.mcmc$freq)
  
})

# Issue #91 remove calls to SETLENGTH by adding GROWABLE vectors
test_that("AMCMC VS MCMC with no adaptation", {
  data(Hald)
  set.seed(42)
  hald.mcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                     data = Hald, method = "MCMC",  burnin.iteration = 200, MCMC.iterations = 0)
  set.seed(42)
  hald.amcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                      data = Hald, method = "AMCMC", burnin.iteration = 200, MCMC.iterations = 0, 
                      GROW = TRUE)
  expect_equal(hald.amcmc$postprobs.MCMC, hald.mcmc$postprobs.MCMC)
  expect_equal(hald.amcmc$postprobs.MCMC, hald.mcmc$postprobs.MCMC)
  expect_equal(hald.amcmc$n.models, hald.mcmc$n.models)
  expect_equal(hald.amcmc$logmarg, hald.mcmc$logmarg)
  expect_equal(hald.amcmc$freq, hald.mcmc$freq)
  
})
