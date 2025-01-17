context("bas.lm Growable vectors")

test_that("Test MCMC with Growable Vectors when not needed", {
  data(Hald)

  bas_hald.old <- bas.lm(Y ~ ., data=Hald, prior="BIC",
                      method="MCMC",  MCMC.init = 10000, burnin = 1000,
                      initprobs = c(1, -.4, .3, 1.0, .8))
  
  bas_hald.grow <- bas.lm(Y ~ ., data=Hald, prior="BIC",
                          method="MCMC_GROWABLE",  MCMC.init = 10000, burnin = 1000,
                          initprobs = c(-.4, .3, 1.5, .8))
  expect_equal(bas_hald.grow$logmarg, bas_hald.old$logmarg)
})