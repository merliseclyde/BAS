#' Coefficients of a Bayesian Model Average object
#'
#' Extract conditional posterior means and standard deviations, marginal
#' posterior means and standard deviations, posterior probabilities, and
#' marginal inclusions probabilities under Bayesian Model Averaging from an
#' object of class 'bas'
#'
#' Calculates posterior means and (approximate) standard deviations of the
#' regression coefficients under Bayesian Model averaging using g-priors and
#' mixtures of g-priors.  Print returns overall summaries. For fully Bayesian
#' methods that place a prior on g, the posterior standard deviations do not
#' take into account full uncertainty regarding g. Will be updated in future
#' releases.
#'
#' @aliases coef.bas coef coefficients coefficients.bas print.coef.bas
#' @param object object of class 'bas' created by BAS
#' @param x object of class 'coef.bas' to print
#' @param n.models Number of top models to report in the printed summary, for
#' coef the default is to use all models.  To extract summaries for the Highest
#' Probability Model, use n.models=1 or estimator="HPM".
#' @param estimator return summaries for a selected model, rather than using
#' BMA.  Options are 'HPM' (highest posterior probability model) ,'MPM' (median
#' probability model), and 'BMA'
#' @param digits number of significant digits to print
#' @param ... other optional arguments
#' @return \code{coefficients} returns an object of class coef.bas with the
#' following:
#'  \item{conditionalmeans}{a matrix with conditional posterior means
#' for each model}
#'  \item{conditionalsd}{ standard deviations for each model }
#'  \item{postmean}{marginal posterior means of each regression coefficient
#' using BMA}
#'  \item{postsd}{marginal posterior standard deviations using BMA}
#'  \item{postne0}{vector of posterior inclusion probabilities, marginal
#' probability that a coefficient is non-zero}
#' @note With highly correlated variables, marginal summaries may not be
#' representative of the joint distribution. Use \code{\link{plot.coef.bas}} to
#' view distributions.  The value reported for the intercept is
#' under the centered parameterization.  Under the  Gaussian error
#' model it will be centered at the sample mean of Y.
#' @author Merlise Clyde \email{clyde@@duke.edu}
#' @seealso \code{\link{bas}}, \code{\link{confint.coef.bas}}
#' @references Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
#' (2005) Mixtures of g-priors for Bayesian Variable Selection.  Journal of the
#' American Statistical Association.  103:410-423.  \cr
#' \doi{10.1198/016214507000001337}
#' @keywords regression
#' @examples
#'
#' data("Hald")
#' hald.gprior =  bas.lm(Y~ ., data=Hald, n.models=2^4, alpha=13,
#'                       prior="ZS-null", initprobs="Uniform", update=10)
#' coef.hald.gprior = coefficients(hald.gprior)
#' coef.hald.gprior
#' plot(coef.hald.gprior)
#' confint(coef.hald.gprior)
#'
#' #Estimation under Median Probability Model
#' coef.hald.gprior = coefficients(hald.gprior, estimator="MPM")
#' coef.hald.gprior
#' plot(coef.hald.gprior)
#' plot(confint(coef.hald.gprior))
#'
#'
#' coef.hald.gprior = coefficients(hald.gprior, estimator="HPM")
#' coef.hald.gprior
#' plot(coef.hald.gprior)
#' confint(coef.hald.gprior)
#'
#' # To add estimation under Best Predictive Model
#'
#'
#' @rdname coef
#' @family bas methods
#' @export
coef.bas <- function(object, n.models, estimator = "BMA", ...) {
  if (estimator == "BPM") {
    stop("Extracting coefficients for the BPM is not implemented yet")
  }
  if (estimator == "MPM") {
    nvar <- object$n.vars - 1
    bestmodel <- (0:nvar)[object$probne0 > .5]
    best <- 1
    models <- rep(0, nvar + 1)
    models[bestmodel + 1] <- 1
    if (sum(models) > 1) {
      
      # fix for issue #39 and #56 
      modelform = as.formula(eval(object$call$formula, parent.frame()))
      environment(modelform) = environment()
      data = eval(object$call$data)
      weights = eval(object$call$weights)
      object <- bas.lm(
        formula=modelform,
        data = data,
        weights = weights,
        n.models = 1,
        alpha = object$g,
        initprobs = object$probne0,
        prior = object$prior,
        modelprior = object$modelprior,
        update = NULL,
        bestmodel = models,
        prob.local = 0.0
      )
    }
  }
  postprobs <- object$postprobs
  if (estimator == "MPM" | estimator == "HPM") {
    n.models <- 1
  }
  if (missing(n.models)) {
    n.models <- length(postprobs)
  }

  topm <- order(-postprobs)[1:n.models]
  postprobs <- postprobs[topm] / sum(postprobs[topm])
  shrinkage <- object$shrinkage[topm]
  conditionalmeans <- list2matrix.bas(object, "mle")[topm, , drop = F]
  conditionalmeans[, -1] <- sweep(conditionalmeans[, -1, drop = F], 1,
    shrinkage,
    FUN = "*"
  )
  postmean <- as.vector(postprobs %*% conditionalmeans)

  # workaround for issue #65
  if (inherits(object, "basglm")) {
     object$prior = object$betaprior$class
  }
  
  conditionalsd <- list2matrix.bas(object, "mle.se")[topm, , drop = F]
  if (!(object$prior == "AIC" | object$prior == "BIC" | object$prior == "IC")) {
    conditionalsd[, -1] <- sweep(conditionalsd[, -1, drop = F], 1,
      sqrt(shrinkage),
      FUN = "*"
    )
  }

  postsd <- sqrt(postprobs %*% conditionalsd^2 +
    postprobs %*% ((sweep(conditionalmeans, 2, postmean, FUN = "-")
    )^2))
  postsd <- as.vector(postsd)
  if (is.null(object$df[topm])) {
    df <- rep(object$n, length(postprobs))
    if (object$prior == "BIC" | object$prior == "AIC" | object$prior == "IC") {
      df <- df - object$rank
    } else {
      df <- df - 1
    }
  } else {
    df <- object$df[topm]
  }

  out <- list(
    postmean = postmean,
    postsd = postsd,
    probne0 = object$probne0,
    conditionalmeans = conditionalmeans,
    conditionalsd = conditionalsd,
    namesx = object$namesx,
    postprobs = postprobs,
    n.vars = object$n.vars,
    n.models = n.models,
    df = df,
    estimator = estimator
  )
  class(out) <- "coef.bas"
  return(out)
}

#' Print coefficients generated from coef.bas
#' @rdname coef
#' @aliases print.coef.bas
#' @family bas coefs
#' @method print coef.bas
#' @export
print.coef.bas <- function(x,
                           digits = max(3, getOption("digits") - 3), ...) {
  out <- cbind(x$postmean, x$postsd, x$probne0)
  dimnames(out) <- list(x$namesx, c("post mean", "post SD", "post p(B != 0)"))

  cat("\n Marginal Posterior Summaries of Coefficients: \n")
  cat("\n Using ", x$estimator, "\n")
  cat("\n Based on the top ", x$n.models, "models \n")
  print.default(format(out, digits = digits),
    print.gap = 2,
    quote = FALSE, ...
  )
  invisible()
}
