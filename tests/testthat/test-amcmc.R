context("AMCMC bas.lm")

test_that("methods", {
  data(Hald)
  set.seed(42)
  hald.amcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                      data = Hald, method = "AMCMC", burnin.iteration = 200, MCMC.iterations = 0)
  set.seed(42)
  hald.mcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                     data = Hald, method = "MCMC",  burnin.iteration = 0, MCMC.iterations = 200)
  expect_equal(hald.amcmc$postprobs.MCMC, hald.mcmc$postprobs.MCMC)
})