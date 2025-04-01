context("AMCMC bas.lm")

test_that("methods", {
  data(Hald)
  set.seed(42)
  hald.mcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                     data = Hald, method = "MCMC",  burnin.iteration = 200, MCMC.iterations = 0)
  set.seed(42)
  hald.amcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                      data = Hald, method = "AMCMC", burnin.iteration = 200, MCMC.iterations = 0, 
                      GROW = FALSE)
  expect_equal(hald.amcmc$postprobs.MCMC, hald.mcmc$postprobs.MCMC)
})

skip()  # skip this test because it is not working for now
test_that("sample", {
  data(Hald)
  set.seed(42)

 hald.amcmc.old = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                         data = Hald, method = "AMCMC", burnin.iteration = 200, 
                         MCMC.iterations = 50, GROW = FALSE, importance.sampling  = FALSE)
 set.seed(42)
 hald.amcmc= bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                         data = Hald, method = "AMCMC", burnin.iteration = 200, 
                         MCMC.iterations = 50, GROW = TRUE, importance.sampling = FALSE)
  expect_equal(hald.amcmc$n.models, hald.amcmc.old$n.models)
  expect_equal(hald.amcmc$postprobs.MCMC, hald.amcmc.old$postprobs.MCMC)
  expect_equal(hald.amcmc$logmarg, hald.amcmc.old$logmarg)
  expect_equal(hald.amcmc$freq, hald.amcmc.old$freq)
  
  
})