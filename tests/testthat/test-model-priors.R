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
