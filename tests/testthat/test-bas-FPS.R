context("bas.lm")

skip_on_os("mac")  # FPS tests are flaky on macOS GH actions
skip_on_arch("aarch64")  # FPS tests are flaky on aarch64 GH actions
# FIXME: Remove the skips above when the underlying issue is resolved
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