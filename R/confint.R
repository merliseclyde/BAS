#' Compute Credible Intervals for BAS regression coefficients from BAS objects
#'
#' Uses Monte Carlo simulations using posterior means and standard deviations
#' of coefficients to generate draws from the posterior distributions and
#' returns highest posterior density (HPD) credible intervals.  If the number
#' of models equals one, then use the t distribution to find intervals.  These
#' currently condition on the estimate of $g$. %% ~~ If necessary, more details
#' than the description above ~~
#'
#' @aliases confint.coef.bas confint
#' @param object a coef.bas object
#' @param parm a specification of which parameters are to be given credible
#' intervals, either a vector of numbers or a vector of names. If missing, all
#' parameters are considered.
#' @param level the probability coverage required
#' @param nsim number of Monte Carlo draws from the posterior distribution.
#' Used when number of models is greater than 1.
#' @param ... other arguments to passed; none currently
#' @return A matrix (or vector) with columns giving lower and upper HPD
#' credible limits for each parameter. These will be labeled as 1-level)/2 and
#' 1 - (1-level)/2 in percent (by default 2.5 and 97.5).
#' @note For mixture of g-priors these are approximate.  This uses Monte Carlo
#' sampling so results may be subject to Monte Carlo variation and larger values
#' of nsim may be needed to reduce variability. %% ~~further notes~~
#' @author Merlise A Clyde
#' @keywords regression
#' @examples
#'
#'
#' data("Hald")
#' hald_gprior <-  bas.lm(Y~ ., data=Hald, alpha=13,
#'                             prior="g-prior")
#' coef_hald <- coef(hald_gprior)
#' confint(coef_hald)
#' confint(coef_hald, approx=FALSE, nsim=5000)
#' # extract just the coefficient of X4
#' confint(coef_hald, parm="X4")
#'
#'
#' @rdname confint.coef
#' @family CI methods
#' @family bas methods
#' @method confint coef.bas
#' @export
confint.coef.bas <- function(object, parm, level = 0.95, nsim = 10000, ...) {
  n.models <- length(object$postprob)
  if (missing(parm)) parm <- 1:object$n.vars
  if (!is.numeric(parm)) parm <- which(object$namesx %in% parm)

  if (n.models > 1) {
    models <- sample(1:n.models, size = nsim, prob = object$postprobs, replace = TRUE)
    means <- object$conditionalmeans[models, parm]
    sd <- object$conditionalsd[models, parm]
    df <- object$df
    if (length(df) == length(object$postprobs)) df <- object$df[models]
    betas <- matrix(rt(nsim * length(subset), df = df),
      nrow = nsim, ncol = length(parm), byrow = FALSE
    )
    betas <- betas * sd + means
    ci <- .HPDinterval(betas, prob = level)
  }
  else {
    df <- sum(object$postprobs * object$df)
    means <- object$postmean[parm]
    sd <- object$postsd[parm]
    tq <- -qt((1 - level) / 2, df = df)
    ci <- cbind(means - tq * sd, means + tq * sd)
  }
  ci <- cbind(ci, object$postmean[parm])
  attr(ci, "Probability") <- level
  attr(ci, "class") <- "confint.bas"
  lower <- paste0(as.character(round(100 * (1 - level) / 2, 4)), "%")
  upper <- paste0(as.character(round(100 * (1 + level) / 2, 4)), "%")
  colnames(ci) <- c(lower, upper, "beta")
  rownames(ci) <- object$namesx[parm]
  return(ci)
}



#' Compute Credible (Bayesian Confidence) Intervals for a BAS predict object
#'
#' Compute credible intervals for in-sample or out of sample prediction or for
#' the regression function
#'
#' This constructs approximate 95 percent Highest Posterior Density intervals
#' for 'pred.bas' objects.  If the estimator is based on model selection, the
#' intervals use a Student t distribution using the estimate of g.  If the
#' estimator is based on BMA, then nsim draws from the mixture of Student t
#' distributions are obtained with the HPD interval obtained from the Monte
#' Carlo draws. %% ~~ If necessary, more details than the description above ~~
#'
#' @param object an object created by \code{\link{predict.bas}}
#' @param parm character variable, "mean" or "pred".  If missing parm='pred'.
#' @param level the nominal level of the (point-wise) credible interval
#' @param nsim number of Monte Carlo simulations for sampling methods with BMA
#' @param ... optional arguments to pass on to next function call; none at this
#' time.
#' @return a matrix with lower and upper level * 100 percent credible intervals
#' for either the mean of the regression function or predicted values.  %%
#' @author Merlise A Clyde
#' @seealso \code{\link{predict.bas}}
#' @keywords regression
#' @examples
#'
#' data("Hald")
#' hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
#' hald.pred = predict(hald.gprior, estimator="BPM", predict=FALSE, se.fit=TRUE)
#' confint(hald.pred, parm="mean")
#' confint(hald.pred)  #default
#' hald.pred = predict(hald.gprior, estimator="BMA", predict=FALSE, se.fit=TRUE)
#' confint(hald.pred)
#'
#'
#' @rdname confint.pred
#' @family bas methods
#' @family CI methods
#' @method confint pred.bas
#' @export
confint.pred.bas <- function(object, parm, level = 0.95, nsim = 10000, ...) {
  if (missing(parm)) parm <- "pred"
  if (parm == "pred") {
    sd <- object$se.pred
  } else {
    sd <- object$se.fit
  }

  if (is.null(sd)) {
    warning("object does not have fitted or prediction standard deviations")
    return()
  }
  if (object$estimator == "BMA") {
    n.models <- length(object$postprob)
    models <- sample(1:n.models, size = nsim, prob = object$postprobs, replace = TRUE)
    means <- object$Ypred[models, ]
    df <- object$df[models]
    #   browser()
    sd <- sd[models, ]
    npred <- length(object$fit)
    pred <- matrix(rt(nsim * npred, df = df),
      nrow = nsim, ncol = npred, byrow = FALSE
    )
    pred <- pred * sd + means
    ci <- .HPDinterval(pred, prob = level)
  }
  else {
    df <- object$df
    means <- object$fit
    tq <- -qt((1 - level) / 2, df = df)
    ci <- cbind(means - tq * sd, means + tq * sd)
  }


  ci <- cbind(ci, object$fit)
  attr(ci, "Probability") <- level
  attr(ci, "class") <- "confint.bas"
  lower <- paste0(as.character(round(100 * (1 - level) / 2, 4)), "%")
  upper <- paste0(as.character(round(100 * (1 + level) / 2, 4)), "%")
  colnames(ci) <- c(lower, upper, parm)
  rownames(ci) <- object$namesx[parm]


  return(ci)
}



#' Plot Bayesian Confidence Intervals
#'
#' Function takes the the output of functions that return credible intervals
#' from BAS objects, and creates a plot of the posterior mean with segments
#' representing the credible interval.  %% ~~ A concise (1-5 lines) description
#' of what the function does. ~~
#'
#' This function takes the HPD intervals or credible intervals created by
#' \code{\link{confint.coef.bas}} or \code{\link{confint.pred.bas}} from BAS
#' objects, and creates a plot of the posterior mean with segments representing
#' the credible interval.  BAS tries to return HPD intervals, and under model
#' averaging these may not be symmetric.  %% ~~ If necessary, more details than
#' the description above ~~
#'
#' @param x the output from \code{\link{confint.coef.bas}} or
#' \code{\link{confint.pred.bas}} containing credible intervals and estimates.
#' @param horizontal orientation of the plot
#' @param ... optional graphical arguments to pass on to plot
#' @return A plot of the credible intervals.
#' @author Merlise A Clyde
#' @seealso \code{\link{confint.coef.bas}}, \code{\link{confint.pred.bas}},
#' \code{\link{coef.bas}}, \code{\link{predict.bas}}, \code{link{bas.lm}}
#' @keywords regression bayesian
#' @examples
#'
#' data(Hald)
#' hald.ZS = bas.lm(Y ~ ., data=Hald, prior="ZS-null", modelprior=uniform())
#' hald.coef = confint(coef(hald.ZS), parm=2:5)
#' plot(hald.coef)
#' plot(hald.coef, horizontal=TRUE)
#' plot(confint(predict(hald.ZS, se.fit=TRUE), parm="mean"))
#'
#' @rdname  plot.confint
#' @method  plot  confint.bas
#' @family bas methods
#' @family CI methods
#' @export
plot.confint.bas <- function(x, horizontal = FALSE, ...) {
  ci <- x #  x is there for generic plot function
  namesx <- rownames(ci)
  y <- ci[, 3]
  x <- 1:nrow(ci)
  xlim <- range(x) + c(-0.5, 0.2)
  ylim <- range(pretty(ci))

  xlab <- "case"

  type <- colnames(ci)[3]
  if (type == "beta") {
    ylab <- bquote(beta)
    xlab <- "coefficient"
  }
  else {
    if (type == "mean") {
      ylab <- bquote(mu)
    } else {
      ylab <- "predicted values"
    }
  }

  not.deg <- ci[, 1] != ci[, 2]
  par(mar = c(4.5, 5, 1, 1), las = 1)
  if (!horizontal) {
    plot(y, pch = 16, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, xaxt = "n", bty = "n", ...)
    axis(1, at = x, labels = namesx, tick = FALSE, ...)
    abline(h = 0, lty = 3, ...)
    arrows(x[not.deg], ci[not.deg, 1], x[not.deg], ci[not.deg, 2], code = 3, angle = 90, length = 0.05, ...)
  }
  ### horizontal layout:
  else {
    plot(
      x = y, y = x, pch = 16, xlim = ylim, ylim = xlim,
      xlab = ylab, ylab = "", yaxt = "n", bty = "n", ...
    )
    axis(2, at = x, labels = namesx, tick = FALSE, ...)
    abline(v = 0, lty = 3, ...)
    arrows(ci[not.deg, 1], x[not.deg], ci[not.deg, 2], x[not.deg],  code = 3, angle = 90, length = 0.05, ...)
  }
  return()
}

.HPDinterval <- function(obj, prob = 0.95, ...) {
  # from library coda but used here so that library
  # does not have to be loaded
  obj <- as.matrix(obj)
  vals <- apply(obj, 2, sort)
  if (!is.matrix(vals)) {
    stop("obj must have nsamp > 1")
  }
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, , drop = FALSE] - vals[init,
    ,
    drop = FALSE
  ], 2, which.min)
  ans <- cbind(vals[cbind(inds, 1:npar)], vals[cbind(inds +
    gap, 1:npar)])
  lower <- as.character(round(100 * (1 - prob) / 2, 2))
  upper <- as.character(round(100 * (prob + 1) / 2, 2))

  dimnames(ans) <- list(colnames(obj), c(lower, upper))
  attr(ans, "Probability") <- gap / nsamp
  ans
}
