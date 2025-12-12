context("bas.lm Growable vectors")

test_that("Test MCMC with Growable Vectors when not needed", {
  data(Hald)


  set.seed(42)
  bas_hald_grow <- bas.lm(Y ~ ., data=Hald, prior="BIC", n.models = 10,
                          method="MCMC",  MCMC.it = 10000, burnin = 1000,
                          GROW = TRUE,
                          initprobs = c(-.4, .3, 1.5, .8))
  
  set.seed(42)
  bas_hald_old <- bas.lm(Y ~ ., data=Hald, prior="BIC", n.models = bas_hald_grow$n.models,
                         method="MCMC",  MCMC.it = 10000, burnin = 1000,
                         GROW = TRUE,
                         initprobs = c(1, -.4, .3, 1.0, .8))
  expect_equal(bas_hald_grow$logmarg, bas_hald_old$logmarg)
  expect_equal(bas_hald_grow$freq, bas_hald_old$freq)
  expect_equal(bas_hald_grow$size, bas_hald_old$size)
  expect_equal(bas_hald_grow$probne0, bas_hald_old$probne0)
  expect_equal(bas_hald_grow$probne0.MCMC, bas_hald_old$probne0.MCMC)
})

test_that("Test BAS with Growable Vectors when not needed", {
  data(Hald)
  
  
  set.seed(42)
  bas_hald_grow <- bas.lm(Y ~ ., data=Hald, prior="BIC",
                          method="BAS", GROW = TRUE,
                          initprobs = c(-.4, .3, 1.5, .8))
  
  set.seed(42)
  bas_hald_old <- bas.lm(Y ~ ., data=Hald, prior="BIC", 
                         method="BAS",  GROW = FALSE,
                         initprobs = c(1, -.4, .3, 1.0, .8))
  expect_equal(bas_hald_grow$n.models, bas_hald_old$n.models)
  expect_equal(bas_hald_grow$logmarg, bas_hald_old$logmarg)
  expect_equal(bas_hald_grow$size, bas_hald_old$size)
  expect_equal(bas_hald_grow$probne0, bas_hald_old$probne0)
  expect_equal(bas_hald_grow$probne0.MCMC, bas_hald_old$probne0.MCMC)
})

#skip() # deprecated
test_that("Test MCMC with Growable Vectors when needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  
  set.seed(42)
  crime.mcmc <-  bas.lm(y ~ ., data=UScrime, n.models=2^16, prior="BIC", 
                        GROW = FALSE,
                        method="MCMC", MCMC.it = 100000, burnin = 1000)
  
  set.seed(42)
  crime.grow = bas.lm(y ~ ., data=UScrime, prior="BIC", 
                      n.models = crime.mcmc$n.models - 100,
                      GROW = TRUE,
                      method="MCMC",MCMC.it = 100000, burnin = 1000)
  
  expect_equal(crime.grow$logmarg, crime.mcmc$logmarg)
  expect_equal(crime.grow$freq, crime.mcmc$freq)
  expect_equal(crime.grow$size, crime.mcmc$size)
  expect_equal(crime.grow$probne0, crime.mcmc$probne0)
})

# passes
test_that("Test AMCMC with Growable Vectors when not needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  
  set.seed(42)
  crime.mcmc <-  bas.lm(y ~ ., data=UScrime, n.models=2^16, prior="BIC", 
                        method="MCMC", GROW = TRUE, MCMC.it = 0, burnin = 100000)
  
  set.seed(42)
  crime.grow = bas.lm(y ~ ., data=UScrime, prior="BIC", 
                      n.models = crime.mcmc$n.models,
                      method="AMCMC",  GROW = TRUE, importance = FALSE,
                      MCMC.it = 0, burnin = 100000)
  
  expect_equal(crime.grow$n.models, crime.mcmc$n.models)
  expect_equal(crime.grow$logmarg, crime.mcmc$logmarg)
  expect_equal(crime.grow$freq, crime.mcmc$freq)
  expect_equal(crime.grow$size, crime.mcmc$size)
  expect_equal(crime.grow$probne0, crime.mcmc$probne0)
})

test_that("Test AMCMC with Growable Vectors when needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  
  set.seed(42)
  crime.mcmc <-  bas.lm(y ~ ., data=UScrime, n.models=2^16, prior="BIC",
                        GROW = TRUE,
                        method="MCMC", MCMC.it = 0, burnin = 100000)
  
  set.seed(42)
  crime.grow <- bas.lm(y ~ ., data=UScrime, prior="BIC", 
                      n.models = crime.mcmc$n.models,
                      method="AMCMC",  GROW = TRUE, importance = FALSE,
                      MCMC.it = 0, burnin = 100000)
  expect_equal(crime.grow$n.models, crime.mcmc$n.models)
  expect_equal(crime.grow$logmarg, crime.mcmc$logmarg)
  expect_equal(crime.grow$freq, crime.mcmc$freq)
  expect_equal(crime.grow$size, crime.mcmc$size)
  expect_equal(crime.grow$probne0, crime.mcmc$probne0)
})

test_that("Test AMCMC with Growable Vectors when not needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  
  set.seed(42)
  crime.mcmc <-  bas.lm(y ~ ., data=UScrime, prior="BIC",
                        n.models=2^16, 
                        method="AMCMC", GROW = FALSE, importance = FALSE,
                        MCMC.it = 0, burnin = 10000)
  
  set.seed(42)
  crime.grow <- bas.lm(y ~ ., data=UScrime, prior="BIC", 
                      n.models = crime.mcmc$n.models,
                      method="AMCMC",  GROW = TRUE, importance = FALSE,
                      MCMC.it = 0, burnin = 10000)
  
  expect_equal(crime.grow$logmarg, crime.mcmc$logmarg)
  expect_equal(crime.grow$freq, crime.mcmc$freq)
  expect_equal(crime.grow$size, crime.mcmc$size)
  expect_equal(crime.grow$probne0, crime.mcmc$probne0)
})

test_that("Test AMCMC with Growable Vectors when needed", {
  # issue #91 implement growable vectors in MCMC_GROWABLE
  data(UScrime, package="MASS")
  UScrime[,-2] <- log(UScrime[,-2])
  
  set.seed(42)
  crime.mcmc <-  bas.lm(y ~ ., data=UScrime, n.models=2^16, prior="BIC",
                        method="AMCMC", GROW = FALSE, MCMC.it = 0, burnin = 100000)
  
  set.seed(42)
  crime.grow <- bas.lm(y ~ ., data=UScrime, prior="BIC", 
                      n.models = crime.mcmc$n.models,
                      method="AMCMC",  GROW = TRUE, importance = FALSE,
                      MCMC.it = 0, burnin = 100000)
  expect_equal(crime.grow$n.models, crime.mcmc$n.models)
  expect_equal(crime.grow$logmarg, crime.mcmc$logmarg)
  expect_equal(crime.grow$freq, crime.mcmc$freq)
  expect_equal(crime.grow$size, crime.mcmc$size)
  expect_equal(crime.grow$probne0, crime.mcmc$probne0)
})
