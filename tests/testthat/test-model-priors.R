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
})
