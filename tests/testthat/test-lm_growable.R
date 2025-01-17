context("bas.lm Growable vectors")

test_that("Test MCMC with Growable Vectors when not needed", {
  data(Hald)


  set.seed(42)
  bas_hald_grow <- bas.lm(Y ~ ., data=Hald, prior="BIC", n.models = 10,
                          method="MCMC_GROWABLE",  MCMC.it = 10000, burnin = 1000,
                          initprobs = c(-.4, .3, 1.5, .8))
  
  set.seed(42)
  bas_hald_old <- bas.lm(Y ~ ., data=Hald, prior="BIC", n.models = bas_hald_grow$n.models,
                         method="MCMC",  MCMC.it = 10000, burnin = 1000,
                         initprobs = c(1, -.4, .3, 1.0, .8))
  expect_equal(bas_hald_grow$logmarg, bas_hald_old$logmarg)
})

test_that("Test MCMC with Growable Vectors when needed", {
  data(Hald)
  expect_no_error(bas.lm(Y ~ ., data=Hald, prior="BIC",
                         method="MCMC", n.models = 2, 
                         MCMC.it = 10000, burnin = 1000,
                         initprobs = c(1, -.4, .3, 1.0, .8)))
  # issue #91 implement growable vectors in MCMC_GROWABLE
  expect_error(bas.lm(Y ~ ., data=Hald, prior="BIC",
                          method="MCMC_GROWABLE", n.models=2,
                          MCMC.it = 10000, burnin = 1000,
                          initprobs = c(-.4, .3, 1.5, .8)))
})
