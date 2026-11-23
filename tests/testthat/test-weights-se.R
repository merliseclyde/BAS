context("standard errors with weights")

# GitHub issue #52: weights (prior weights in linear models and the IRLS
# weights in GLMs) were not accounted for in the covariance matrix used for
# standard errors.  With the BIC prior there is no shrinkage and
# the posterior means and standard deviations given a model agree with the
# MLEs and their standard errors, so that results may be compared directly
# to lm() and glm().

sim_weights_data <- function() {
  set.seed(52)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$w <- sample(1:4, n, replace = TRUE) + runif(n)   # non-integer weights
  d$wi <- sample(1:4, n, replace = TRUE)             # integer weights
  d$y <- 1 + d$x1 - .5 * d$x2 + rnorm(n) / sqrt(d$w)
  d$z <- rbinom(n, 1, plogis(.3 + d$x1))
  d$cnt <- rpois(n, exp(.3 + .5 * d$x1))
  d
}

newdata_weights <- data.frame(
  x1 = c(-1, 0, 2),
  x2 = c(.5, 0, -1),
  x3 = c(0, 1, 3)
)

# formula for the model with the variables in `vars`
model_formula <- function(vars, response) {
  if (length(vars) == 0) {
    as.formula(paste(response, "~ 1"))
  } else {
    reformulate(vars, response)
  }
}

model_vars <- function(object, model) {
  object$namesx[model[-1] + 1]
}

test_that("weighted regression: coefficients and SEs agree with lm", {
  d <- sim_weights_data()
  # bas.lm centers the predictors at their weighted means
  dc <- d
  for (v in c("x1", "x2", "x3")) {
    dc[[v]] <- d[[v]] - weighted.mean(d[[v]], d$w)
  }

  for (method in c("BAS", "deterministic", "MCMC")) {
    for (pivot in c(TRUE, FALSE)) {
      wt_bas <- bas.lm(y ~ x1 + x2 + x3,
        data = d, weights = w,
        prior = "BIC", modelprior = uniform(),
        method = method, pivot = pivot
      )
      expect_equal(wt_bas$weights, d$w)

      for (m in seq_along(wt_bas$which)) {
        vars <- model_vars(wt_bas, wt_bas$which[[m]])
        wt_lm <- lm(model_formula(vars, "y"), data = dc, weights = w)
        expect_equal(wt_bas$mle[[m]], unname(coef(wt_lm)))
        expect_equal(wt_bas$mle.se[[m]], unname(sqrt(diag(vcov(wt_lm)))))

        # slopes and SEs do not depend on centering
        wt_lm_raw <- lm(model_formula(vars, "y"), data = d, weights = w)
        expect_equal(wt_bas$mle[[m]][-1], unname(coef(wt_lm_raw))[-1])
        expect_equal(
          wt_bas$mle.se[[m]][-1],
          unname(sqrt(diag(vcov(wt_lm_raw))))[-1]
        )
      }
    }
  }
})

test_that("weighted regression: se.fit and se.pred agree with predict.lm", {
  d <- sim_weights_data()
  wt_bas <- bas.lm(y ~ x1 + x2 + x3,
    data = d, weights = w,
    prior = "BIC", modelprior = uniform()
  )

  for (estimator in c("HPM", "MPM", "BPM")) {
    # fitted values
    wt_pred <- predict(wt_bas, estimator = estimator, se.fit = TRUE)
    wt_lm <- lm(model_formula(model_vars(wt_bas, wt_pred$bestmodel), "y"),
      data = d, weights = w
    )
    lm_pred <- predict(wt_lm, se.fit = TRUE)
    expect_equal(as.vector(wt_pred$fit), unname(lm_pred$fit))
    expect_equal(as.vector(wt_pred$se.fit), unname(lm_pred$se.fit))

    # new data
    wt_pred <- predict(wt_bas,
      newdata = newdata_weights,
      estimator = estimator, se.fit = TRUE
    )
    wt_lm <- lm(model_formula(model_vars(wt_bas, wt_pred$bestmodel), "y"),
      data = d, weights = w
    )
    lm_pred <- predict(wt_lm, newdata = newdata_weights, se.fit = TRUE)
    expect_equal(as.vector(wt_pred$fit), unname(lm_pred$fit))
    expect_equal(as.vector(wt_pred$se.fit), unname(lm_pred$se.fit))
    # prediction of a new observation with unit weight
    expect_equal(
      as.vector(wt_pred$se.pred),
      unname(sqrt(lm_pred$se.fit^2 + lm_pred$residual.scale^2))
    )
    # credible intervals agree with confidence and prediction intervals
    lm_ci <- predict(wt_lm, newdata = newdata_weights, interval = "confidence")
    expect_equivalent(
      unclass(confint(wt_pred, parm = "mean"))[, 1:2],
      lm_ci[, c("lwr", "upr")]
    )
    # predict.lm warns that it assumes unit weights for new observations
    lm_pi <- suppressWarnings(
      predict(wt_lm, newdata = newdata_weights, interval = "prediction")
    )
    expect_equivalent(
      unclass(confint(wt_pred, parm = "pred"))[, 1:2],
      lm_pi[, c("lwr", "upr")]
    )
  }

  # model specific standard errors under BMA
  for (newdata in list(NULL, newdata_weights)) {
    wt_pred <- predict(wt_bas,
      newdata = newdata, estimator = "BMA",
      se.fit = TRUE
    )
    lm_newdata <- if (is.null(newdata)) d else newdata
    for (m in seq_along(wt_pred$best)) {
      vars <- model_vars(wt_bas, wt_bas$which[[wt_pred$best[m]]])
      wt_lm <- lm(model_formula(vars, "y"), data = d, weights = w)
      lm_pred <- predict(wt_lm, newdata = lm_newdata, se.fit = TRUE)
      expect_equal(wt_pred$se.fit[m, ], unname(lm_pred$se.fit))
      expect_equal(
        wt_pred$se.pred[m, ],
        unname(sqrt(lm_pred$se.fit^2 + lm_pred$residual.scale^2))
      )
    }
  }
})

test_that("weights equal to one give the same SEs as no weights", {
  d <- sim_weights_data()
  d$one <- 1
  no_wt <- bas.lm(y ~ x1 + x2 + x3, data = d, prior = "ZS-null")
  wt_one <- bas.lm(y ~ x1 + x2 + x3, data = d, weights = one, prior = "ZS-null")
  for (estimator in c("HPM", "MPM", "BMA")) {
    p0 <- predict(no_wt,
      newdata = newdata_weights, estimator = estimator,
      se.fit = TRUE
    )
    p1 <- predict(wt_one,
      newdata = newdata_weights, estimator = estimator,
      se.fit = TRUE
    )
    expect_equal(p0$se.fit, p1$se.fit)
    expect_equal(p0$se.pred, p1$se.pred)
    expect_equal(p0$se.bma.fit, p1$se.bma.fit)
  }
})

test_that("weighted regression with shrinkage: se.fit under the g-prior", {
  d <- sim_weights_data()
  n <- nrow(d)
  g <- as.numeric(n)
  wt_bas <- bas.lm(y ~ x1 + x2 + x3,
    data = d, weights = w,
    prior = "g-prior", alpha = g, modelprior = uniform()
  )
  wt_pred <- predict(wt_bas,
    newdata = newdata_weights, estimator = "HPM",
    se.fit = TRUE
  )
  vars <- model_vars(wt_bas, wt_pred$bestmodel)
  wt_lm <- lm(model_formula(vars, "y"), data = d, weights = w)

  # posterior under the g-prior with predictors centered at weighted means:
  # Var(intercept) = sigma^2/sum(w)
  # Cov(beta) = g/(1 + g) sigma^2 (Xc^T W Xc)^{-1}
  shrinkage <- g / (1 + g)
  expect_equal(wt_bas$shrinkage[wt_pred$best], shrinkage)

  xbar <- sapply(d[, vars, drop = FALSE], weighted.mean, w = d$w)
  Xc <- sweep(as.matrix(d[, vars, drop = FALSE]), 2, xbar)
  Xnew <- sweep(as.matrix(newdata_weights[, vars, drop = FALSE]), 2, xbar)
  V <- solve(crossprod(sqrt(d$w) * Xc))
  scale <- 1 / sum(d$w) + shrinkage * rowSums((Xnew %*% V) * Xnew)

  ssy <- sum(d$w * (d$y - weighted.mean(d$y, d$w))^2)
  R2 <- summary(wt_lm)$r.squared
  sigma2 <- ssy * (1 - shrinkage * R2) / (n - 1)

  expect_equal(as.vector(wt_pred$se.fit), unname(sqrt(sigma2 * scale)))
  expect_equal(as.vector(wt_pred$se.pred), unname(sqrt(sigma2 * (1 + scale))))
})

test_that("logistic and Poisson regression: SEs agree with glm", {
  d <- sim_weights_data()

  # Tighter convergence than the default glm.control() (epsilon = 1e-8):
  # predict.glm uses the IRLS weights from the last iteration, evaluated one
  # iterate before the final estimate, while predict.basglm recomputes them
  # at the MLE, so at the default the two differ by up to about 1e-4;
  # with epsilon = 1e-12 the largest relative difference is about 5e-7
  ctl <- glm.control(epsilon = 1e-12, maxit = 100)

  fits <- list(
    logistic = list(
      bas = bas.glm(z ~ x1 + x2 + x3,
        data = d, family = binomial(),
        betaprior = bic.prior(), modelprior = uniform(),
        control = ctl
      ),
      glm = function(f) glm(f, data = d, family = binomial(),
        control = ctl
      ),
      response = "z"
    ),
    weighted_logistic = list(
      bas = bas.glm(z ~ x1 + x2 + x3,
        data = d, weights = wi, family = binomial(),
        betaprior = bic.prior(), modelprior = uniform(),
        control = ctl
      ),
      glm = function(f) glm(f, data = d, weights = wi, family = binomial(),
        control = ctl
      ),
      response = "z"
    ),
    poisson = list(
      bas = bas.glm(cnt ~ x1 + x2 + x3,
        data = d, family = poisson(),
        betaprior = bic.prior(), modelprior = uniform(),
        control = ctl
      ),
      glm = function(f) glm(f, data = d, family = poisson(),
        control = ctl
      ),
      response = "cnt"
    ),
    weighted_poisson = list(
      bas = bas.glm(cnt ~ x1 + x2 + x3,
        data = d, weights = wi, family = poisson(),
        betaprior = bic.prior(), modelprior = uniform(),
        control = ctl
      ),
      glm = function(f) glm(f, data = d, weights = wi, family = poisson(),
        control = ctl
      ),
      response = "cnt"
    )
  )

  tol <- 1e-6

  for (fit in fits) {
    glm_bas <- fit$bas
    response <- fit$response

    # coefficients and standard errors for every model
    for (m in seq_along(glm_bas$which)) {
      vars <- model_vars(glm_bas, glm_bas$which[[m]])
      glm_m <- fit$glm(model_formula(vars, response))
      expect_equal(glm_bas$mle[[m]], unname(coef(glm_m)), tolerance = tol)
      expect_equal(glm_bas$mle.se[[m]], unname(sqrt(diag(vcov(glm_m)))),
        tolerance = tol
      )
    }

    for (type in c("link", "response")) {
      for (estimator in c("HPM", "BPM")) {
        # fitted values
        glm_pred <- predict(glm_bas,
          estimator = estimator, type = type,
          se.fit = TRUE
        )
        vars <- model_vars(glm_bas, glm_pred$bestmodel)
        ref <- predict(fit$glm(model_formula(vars, response)),
          type = type, se.fit = TRUE
        )
        expect_equal(as.vector(glm_pred$fit), unname(ref$fit), tolerance = tol)
        expect_equal(as.vector(glm_pred$se.fit), unname(ref$se.fit),
          tolerance = tol
        )

        # new data
        glm_pred <- predict(glm_bas,
          newdata = newdata_weights,
          estimator = estimator, type = type, se.fit = TRUE
        )
        vars <- model_vars(glm_bas, glm_pred$bestmodel)
        ref <- predict(fit$glm(model_formula(vars, response)),
          newdata = newdata_weights, type = type, se.fit = TRUE
        )
        expect_equal(as.vector(glm_pred$fit), unname(ref$fit), tolerance = tol)
        expect_equal(as.vector(glm_pred$se.fit), unname(ref$se.fit),
          tolerance = tol
        )
      }

      # model specific standard errors under BMA
      glm_pred <- predict(glm_bas,
        newdata = newdata_weights,
        estimator = "BMA", type = type, se.fit = TRUE
      )
      mu <- matrix(NA, length(glm_pred$best), nrow(newdata_weights))
      for (m in seq_along(glm_pred$best)) {
        vars <- model_vars(glm_bas, glm_bas$which[[glm_pred$best[m]]])
        ref <- predict(fit$glm(model_formula(vars, response)),
          newdata = newdata_weights, type = type, se.fit = TRUE
        )
        mu[m, ] <- ref$fit
        expect_equal(glm_pred$se.fit[m, ], unname(ref$se.fit), tolerance = tol)
      }
      # total variance = E[Var | model] + Var(E[ | model])
      expect_equal(as.vector(glm_pred$fit),
        as.vector(glm_pred$postprobs %*% mu),
        tolerance = tol
      )
      se_bma <- sqrt(
        as.vector(glm_pred$postprobs %*% glm_pred$se.fit^2) +
          as.vector(glm_pred$postprobs %*%
            sweep(mu, 2, as.vector(glm_pred$fit))^2)
      )
      expect_equal(glm_pred$se.bma.fit, se_bma, tolerance = tol)
    }
  }
})

test_that("weights that are a variable in the data frame work with MPM", {
  d <- sim_weights_data()
  wt_bas <- bas.lm(y ~ x1 + x2 + x3,
    data = d, weights = w,
    prior = "BIC", modelprior = uniform()
  )
  expect_no_error(predict(wt_bas, estimator = "MPM", se.fit = TRUE))
  expect_no_error(coef(wt_bas, estimator = "MPM"))
})

test_that("MPM that is the null model returns the weighted mean", {
  set.seed(3)
  n <- 40
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), w = runif(n, 1, 3))
  d$y <- 2 + rnorm(n) / sqrt(d$w)
  null_bas <- bas.lm(y ~ x1 + x2,
    data = d, weights = w,
    prior = "BIC", modelprior = uniform()
  )
  skip_if(any(null_bas$probne0[-1] > .5))
  null_pred <- predict(null_bas, estimator = "MPM", se.fit = TRUE)
  lm_pred <- predict(lm(y ~ 1, data = d, weights = w), se.fit = TRUE)
  expect_equal(as.vector(null_pred$fit), unname(lm_pred$fit))
  expect_equal(as.vector(null_pred$se.fit), unname(lm_pred$se.fit))
})
