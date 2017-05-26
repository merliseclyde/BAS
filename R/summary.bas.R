#' Print a Summary of Bayesian Model Averaging objects from BAS
#'
#' \code{summary} and \code{print} methods for Bayesian model averaging objects
#' created by \code{bas} Bayesian Adaptive Sampling
#'
#' The print methods display a view similar to \code{print.lm} .  The summary
#' methods display a view specific to Bayesian model averaging giving the top 5
#' highest probability models represented by their inclusion indicators.
#' Summaries of the models include the Bayes Factor (BF) of each model to the
#' model with the largest marginal likelihood, the posterior probabilty of the
#' models, R2, dim (which includes the intercept) and the log of the marginal
#' likelihood.
#'
#' @aliases print.bas print
#' @param x object of class 'bas'
#' @param digits optional number specifying the number of digits to display
#' @param ... other parameters to be passed to \code{print.default}
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{coef.bas}}
#' @keywords print regression
#' @examples
#'
#' library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
#' print(crime.bic)
#' summary(crime.bic)
#' @rdname print
#' @method print bas
#' @export
#'
print.bas = function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Marginal Posterior Inclusion Probabilities: \n")
    out = x$probne0
    names(out) = x$namesx
    print.default(format(out, digits = digits), print.gap = 2,
        quote = FALSE, ...)
    invisible()
}



#' Summaries of Bayesian Model Averaging objects from BAS
#'
#' \code{summary} and \code{print} methods for Bayesian model averaging objects
#' created by \code{bas} Bayesian Adaptive Sampling
#'
#' The print methods display a view similar to \code{print.lm} .  The summary
#' methods display a view specific to Bayesian model averaging giving the top 5
#' highest probability models represented by their inclusion indicators.
#' Summaries of the models include the Bayes Facotr (BF) of each model to the
#' model with the largest marginal likelihood, the posterior probabilty of the
#' models, R2, dim (which includes the intercept) and the log of the marginal
#' likelihood.
#'
#' @aliases summary.bas summary
#' @param object object of class 'bas'
#' @param n.models optional number specifying the number of best models to
#' display in summary
#' @param ... other parameters to be passed to \code{summary.default}
#' @author Merlise Clyde \email{clyde@@duke.edu}
#' @seealso \code{\link{coef.bas}}
#' @keywords print regression
#' @examples
#'
#' library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
#' print(crime.bic)
#' summary(crime.bic)
#' @rdname summary
#' @method summary bas
#' @family bas methods
#' @export
summary.bas = function(object, n.models = 5, ...) {
  best = order(-object$postprobs)
  n.models = min(n.models, length(best))
  best = best[1:n.models]
  x = cbind(list2matrix.which(object, best),
               exp(object$logmarg[best] - max(object$logmarg[best])),
               round(object$postprobs[best], 4),
               round(object$R2[best], 4),
               object$size[best],
               object$logmarg[best])
  x = t(x)
  x = cbind(NA, x)

  x[1:object$n.vars, 1] = object$probne0
  colnames(x) = c("P(B != 0 | Y)", paste("model", 1:n.models))
  rownames(x) = c(object$namesx, "BF", "PostProbs", "R2", "dim", "logmarg")
  return(x)
  }


