#' Update BAS object using a new prior
#'
#' Update a BMA object using a new prior distribution on the coefficients.
#'
#' Recomputes the marginal likelihoods for the new methods for models already
#' sampled in current object.
#'
#' @aliases update update.bas
#' @param object BMA object to update
#' @param newprior Update posterior model probabilities, probne0, shrinkage,
#' logmarg, etc, using prior based on newprior.  See \code{\link{bas}} for
#' available methods
#' @param alpha optional new value of hyperparameter in prior for method
#' @param ... optional arguments
#' @return A new object of class BMA
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{bas}} for available methods and choices of alpha
#' @references Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive
#' Sampling for Variable Selection and Model Averaging. Journal of
#' Computational Graphics and Statistics.  20:80-101 \cr
#' \url{https://dx.doi.org/10.1198/jcgs.2010.09049}
#' @keywords regression
#' @examples
#'
#' \dontrun{
#' library(MASS)
#' data(UScrime)
#' UScrime[,-2] <- log(UScrime[,-2])
#' crime.bic <-  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
#' crime.ebg <- update(crime.bic, newprior="EB-global")
#' crime.zs <- update(crime.bic, newprior="ZS-null")
#' }
#'
#' @rdname update
#' @method update bas
#' @family bas methods
#' @export
update.bas <- function(object, newprior, alpha = NULL, ...) {
  method.num <- switch(newprior,
    "g-prior" = 0,
    "hyper-g" = 1,
    "EB-local" = 2,
    "BIC" = 3,
    "ZS-null" = 4,
    "ZS-full" = 5,
    "hyper-g-laplace" = 6,
    "AIC" = 7,
    "EB-global" = 2,
    "hyper-g-n" = 8,
    "JZS" = 9,
  )
  if (is.null(alpha) &&
    (method.num == 0 || method.num == 1 || method.num == 6 || method.num == 8)) {
    stop(paste("Must specify a value of alpha for", newprior))
  }

  if (is.null(alpha)) alpha <- 0.0
  object$alpha <- alpha

  if (newprior == "EB-global") {
    object <- EB.global(object)
  } else {
    object$prior <- newprior
    SSY <- sum((object$Y - mean(object$Y))^2)
    R2Full <- summary(lm(object$Y ~ object$X[, -1]))$r.squared
    logmarg <- object$logmarg
    shrinkage <- object$shrinkage
    tmp <- .C(C_gexpectations_vect,
      nmodels = as.integer(length(object$which)),
      p = as.integer(object$n.vars), pmodel = as.integer(object$rank),
      nobs = as.integer(object$n), R2 = object$R2, alpha = as.double(alpha),
      method = as.integer(method.num), RSquareFull = as.double(R2Full), SSY = as.double(SSY),
      logmarg = logmarg, shrinkage = shrinkage
    )
    object$logmarg <- tmp$logmarg
    object$shrinkage <- tmp$shrinkage
    object$postprobs <- exp(object$logmarg - min(object$logmarg)) * object$priorprobs
    object$postprobs <- object$postprobs / sum(object$postprobs)
    which <- which.matrix(object$which, object$n.vars)
    object$probne0 <- object$postprobs %*% which
  }
  return(object)
}
