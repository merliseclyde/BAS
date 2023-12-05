context("model priors")

test_that("bernoulli", {
  data(Hald)
  hald_unif <- bas.lm(Y ~ .,
    data = Hald, prior = "g-prior",
    modelprior = uniform()
  )
  hald_ber <- bas.lm(Y ~ .,
    data = Hald, prior = "g-prior",
    modelprior = Bernoulli(probs = .5)
  )
  expect_equal(hald_unif$probne0, hald_ber$probne0)
  hald_ber1 <- bas.lm(Y ~ .,
                     data = Hald, prior = "g-prior",
                     modelprior = Bernoulli(probs = .4)
  )
  hald_ber2 <- bas.lm(Y ~ X1 + X2 + X3 + X4,
                      data = Hald, prior = "g-prior",
                      modelprior = Bernoulli(probs = rep(.4, 4))
  )
  expect_equal(hald_ber1$probne0, hald_ber2$probne0)
})

test_that("truncated prior", {
hald_tr_power <- bas.lm(Y ~ .,
                    data = Hald, prior = "g-prior",
                    modelprior = tr.power.prior(kappa=2, 2))
expect_equal(1, sum(hald_tr_power$postprobs))
expect_no_error(expect_equal(0, sum(hald_tr_power$postprobs <= 0.0)))

})

test_that("Bernoulli hereditary prior", {
  expect_error(bas.lm(Y ~ .,
                      data = Hald, prior = "g-prior",
                      modelprior = Bernoulli.heredity(.5, NULL))
  )
})

test_that("Always include changes model-prior", {
  
  res <- bas.lm(Y ~ ., data = Hald, modelprior = beta.binomial(1, 1), 
                       include.always = ~ 1 + X1 + X2)
  ord <- order(lengths(res$which))
  expect_equal(
    object = res$priorprobs[ord],
    expected = c(0.333333333333333, 0.166666666666667, 0.166666666666667, 0.333333333333333)
  )
})

test_that("Check if model-prior is of class prior", {
  
  expect_error(bas.lm(Y ~ ., data = Hald, modelprior = c(.25, .25, .25, .25), 
                include.always = ~ 1 + X1 + X2))
  
  expect_error(bas.glm(Y ~ ., data = Hald, modelprior = c(.25, .25, .25, .25), 
                      include.always = ~ 1 + X1 + X2, family=poisson()))
})
