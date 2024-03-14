context("AMCMC bas.lm")

test_that("methods", {
  data(Hald)
  set.seed(42)
  hald.amcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                      data = Hald, method = "AMCMC", MCMC.iterations = 200)
  set.seed(42)
  hald.mcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                     data = Hald, method = "MCMC", MCMC.iterations = 400)
  expect_equal(hald.amcmc$probs.MCMC, hald.mcmc$probs.MCMC)
})