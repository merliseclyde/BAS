context("fitting a single model with n.models = 1")

# GitHub issue #107: with n.models = 1 and bestmodel, the MCMC methods should
# return results for the single model bestmodel, as with method = "BAS".
# A single model is now fit with method = "BAS" for every method, as repeated
# MCMC sampling at one model adds no information.  Previously the growable
# samplers wrote past the end of their output vectors and the samplers
# without GROW stopped before any iterations.

sim_single_model_data <- function() {
  set.seed(42)
  n <- 120
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- rbinom(n, 1, plogis(1.2 * d$x1))
  d$yc <- 1 + d$x1 + rnorm(n)
  d
}

expect_same_single_model <- function(fit, ref) {
  expect_length(fit$which, 1)
  expect_equal(fit$n.models, 1)
  expect_equal(fit$which, ref$which)
  expect_equal(fit$mle, ref$mle)
  expect_equal(fit$mle.se, ref$mle.se)
  expect_equal(fit$logmarg, ref$logmarg)
  expect_equal(fit$postprobs, ref$postprobs)
  expect_equal(fit$probne0, ref$probne0)
  fit_coef <- coef(fit)
  ref_coef <- coef(ref)
  expect_equal(fit_coef$postmean, ref_coef$postmean)
  expect_equal(fit_coef$postsd, ref_coef$postsd)
}

# exact example from the issue
test_that("bas.glm with method = 'MCMC' and n.models = 1 (issue #107)", {
  set.seed(42)
  n <- 120
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- rbinom(n, 1, plogis(1.2 * d$x1))

  # a single known model to evaluate: intercept + x1
  best <- c(1, 1, 0, 0)

  # works with method = "BAS"
  fit_bas <- bas.glm(y ~ x1 + x2 + x3,
    data = d, family = binomial(),
    method = "BAS", n.models = 1, bestmodel = best
  )
  expect_equal(fit_bas$which, list(c(0, 1)))

  # failed with method = "MCMC" and n.models = 1
  expect_no_error(
    fit_mcmc <- bas.glm(y ~ x1 + x2 + x3,
      data = d, family = binomial(),
      method = "MCMC", n.models = 1, bestmodel = best
    )
  )
  expect_same_single_model(fit_mcmc, fit_bas)
})

test_that("bas.glm: every sampling method fits a single model", {
  d <- sim_single_model_data()
  best <- c(1, 1, 0, 0)
  ref <- bas.glm(y ~ x1 + x2 + x3,
    data = d, family = binomial(),
    method = "BAS", n.models = 1, bestmodel = best
  )
  for (method in c("MCMC", "MCMC+BAS")) {
    for (GROW in c(TRUE, FALSE)) {
      fit <- bas.glm(y ~ x1 + x2 + x3,
        data = d, family = binomial(),
        method = method, n.models = 1, bestmodel = best, GROW = GROW
      )
      expect_same_single_model(fit, ref)
    }
  }
  # the single model is fit with method = "BAS"
  fit <- bas.glm(y ~ x1 + x2 + x3,
    data = d, family = binomial(),
    method = "MCMC", n.models = 1, bestmodel = best
  )
  expect_equal(fit$call$method, "BAS")
})

test_that("bas.lm: every sampling method fits a single model", {
  d <- sim_single_model_data()
  best <- c(1, 1, 0, 0)
  ref <- bas.lm(yc ~ x1 + x2 + x3,
    data = d,
    method = "BAS", n.models = 1, bestmodel = best
  )
  for (method in c("MCMC", "MCMC+BAS", "AMCMC")) {
    fit <- bas.lm(yc ~ x1 + x2 + x3,
      data = d,
      method = method, n.models = 1, bestmodel = best
    )
    expect_same_single_model(fit, ref)
  }
  fit <- bas.lm(yc ~ x1 + x2 + x3,
    data = d,
    method = "AMCMC", importance.sampling = TRUE,
    n.models = 1, bestmodel = best
  )
  expect_same_single_model(fit, ref)
  for (method in c("BAS", "AMCMC")) {
    fit <- bas.lm(yc ~ x1 + x2 + x3,
      data = d,
      method = method, GROW = FALSE, n.models = 1, bestmodel = best
    )
    expect_same_single_model(fit, ref)
  }
})

test_that("single model with include.always and force.heredity", {
  d <- sim_single_model_data()
  # columns: (Intercept), x1, x2, x3, x1:x2
  best <- c(1, 1, 1, 1, 1)
  ref <- bas.lm(yc ~ x1 * x2 + x3,
    data = d, method = "BAS", n.models = 1, bestmodel = best,
    include.always = ~ x3, force.heredity = TRUE
  )
  expect_equal(ref$which, list(0:4))
  for (method in c("MCMC", "MCMC+BAS", "AMCMC")) {
    fit <- bas.lm(yc ~ x1 * x2 + x3,
      data = d, method = method, n.models = 1, bestmodel = best,
      include.always = ~ x3, force.heredity = TRUE
    )
    expect_same_single_model(fit, ref)
  }
  ref <- bas.glm(y ~ x1 * x2 + x3,
    data = d, family = binomial(),
    method = "BAS", n.models = 1, bestmodel = best,
    include.always = ~ x3, force.heredity = TRUE
  )
  expect_equal(ref$which, list(0:4))
  for (method in c("MCMC", "MCMC+BAS")) {
    fit <- bas.glm(y ~ x1 * x2 + x3,
      data = d, family = binomial(),
      method = method, n.models = 1, bestmodel = best,
      include.always = ~ x3, force.heredity = TRUE
    )
    expect_same_single_model(fit, ref)
  }
})

test_that("MCMC with a small n.models > 1 grows the output", {
  # n.models = 2 or 3 is the initial size of the output and wrote past the end
  # (attempt to set index 2/2 in SET_VECTOR_ELT)
  d <- sim_single_model_data()
  best <- c(1, 1, 0, 0)
  for (k in 2:3) {
    for (method in c("MCMC", "AMCMC")) {
      set.seed(1)
      fit_k <- bas.lm(yc ~ x1 + x2 + x3,
        data = d,
        method = method, n.models = k, bestmodel = best
      )
      set.seed(1)
      fit <- bas.lm(yc ~ x1 + x2 + x3,
        data = d,
        method = method, bestmodel = best
      )
      expect_equal(fit_k$n.models, fit$n.models)
      expect_equal(fit_k$probne0, fit$probne0)
    }
    set.seed(1)
    fit_k <- bas.glm(y ~ x1 + x2 + x3,
      data = d, family = binomial(),
      method = "MCMC", n.models = k, bestmodel = best
    )
    set.seed(1)
    fit <- bas.glm(y ~ x1 + x2 + x3,
      data = d, family = binomial(),
      method = "MCMC", bestmodel = best
    )
    expect_equal(fit_k$n.models, fit$n.models)
    expect_equal(fit_k$probne0, fit$probne0)
  }
})

test_that("output vectors grow from a small initial size", {
  # (int) (expand * n) did not increase small n (e.g. 1.05 * 2 = 2)
  d <- sim_single_model_data()
  set.seed(1)
  fit_small <- bas.lm(yc ~ x1 + x2 + x3,
    data = d, method = "MCMC",
    n.models.init = 1, MCMC.iterations = 2000
  )
  set.seed(1)
  fit <- bas.lm(yc ~ x1 + x2 + x3,
    data = d, method = "MCMC",
    MCMC.iterations = 2000
  )
  expect_equal(fit_small$n.models, fit$n.models)
  expect_equal(fit_small$probne0, fit$probne0)

  set.seed(1)
  fit_small <- bas.glm(y ~ x1 + x2 + x3,
    data = d, family = binomial(), method = "MCMC",
    n.models.init = 1, MCMC.iterations = 2000
  )
  set.seed(1)
  fit <- bas.glm(y ~ x1 + x2 + x3,
    data = d, family = binomial(), method = "MCMC",
    MCMC.iterations = 2000
  )
  expect_equal(fit_small$n.models, fit$n.models)
  expect_equal(fit_small$probne0, fit$probne0)
})

test_that("growth with expand = 1 and up to the bound on the number of models", {
  # expand = 1 grows by one slot; a large expand is capped at
  # min(2^p, burnin.iterations + MCMC.iterations + 1)
  d <- sim_single_model_data()
  for (expand in c(1, 10)) {
    for (method in c("MCMC", "AMCMC")) {
      set.seed(1)
      fit_grow <- bas.lm(yc ~ x1 + x2 + x3,
        data = d, method = method,
        n.models.init = 2, expand = expand
      )
      set.seed(1)
      fit <- bas.lm(yc ~ x1 + x2 + x3, data = d, method = method)
      expect_equal(fit_grow$n.models, fit$n.models)
      expect_equal(fit_grow$probne0, fit$probne0)
    }
    set.seed(1)
    fit_grow <- bas.glm(y ~ x1 + x2 + x3,
      data = d, family = binomial(), method = "MCMC",
      n.models.init = 2, expand = expand
    )
    set.seed(1)
    fit <- bas.glm(y ~ x1 + x2 + x3,
      data = d, family = binomial(), method = "MCMC"
    )
    expect_equal(fit_grow$n.models, fit$n.models)
    expect_equal(fit_grow$probne0, fit$probne0)
  }
})

test_that("AMCMC grows the output after the burn-in", {
  # with a short burn-in new models are found in the second stage of AMCMC
  d <- sim_single_model_data()
  set.seed(1)
  fit_grow <- bas.lm(yc ~ x1 + x2 + x3,
    data = d, method = "AMCMC",
    n.models.init = 2, burnin.iterations = 1, MCMC.iterations = 2000
  )
  set.seed(1)
  fit <- bas.lm(yc ~ x1 + x2 + x3,
    data = d, method = "AMCMC",
    burnin.iterations = 1, MCMC.iterations = 2000
  )
  expect_gt(fit_grow$n.models, 2)
  expect_equal(fit_grow$n.models, fit$n.models)
  expect_equal(fit_grow$probne0, fit$probne0)
})
