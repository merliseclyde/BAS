context("bas.lm")

skip()  # FPS tests are flaky due to difference in calculations of small probabilities on different systems

# FIXME: Remove the skips above when the underlying issue is resolved due to normalization of very small probabilities
test_that("FPS enumerate", {
  data("Hald")

  set.seed(42)
  hald.bas = bas.lm(Y ~ .,
         prior = "BIC", method = "MCMC+BAS", 
         burnin.iterations = 1000, n.models=2^4,
         modelprior = uniform(), data = Hald,
         FPS = "none")

  set.seed(42)
  hald.bas.bayes = bas.lm(Y ~ .,
                  prior = "BIC", method = "MCMC+BAS", 
                  burnin.iterations = 1000, n.models=2^4,
                  modelprior = uniform(), data = Hald,
                  FPS = "Bayes_HT")

  expect_equal(hald.bas$postprobs, hald.bas.bayes$postprobs)
})