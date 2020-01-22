#' Find the global Empirical Bayes estimates for BMA
#'
#' Finds the global Empirical Bayes estimates of g in Zellner's g-prior and
#' model probabilities
#'
#' Uses the EM algorithm in Liang et al to estimate the type II MLE of g in
#' Zellner's g prior
#'
#' @aliases EB.global EB.global.bas
#' @param object A 'bas' object created by \code{\link{bas}}
#' @param tol tolerance for estimating g
#' @param g.0 initial value for g
#' @param max.iterations Maximum number of iterations for the EM algorithm
#' @return An object of class 'bas' using Zellner's g prior with an estimate of
#' g based on all models
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{bas}}, \code{\link{update}}
#' @references Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
#' (2008) Mixtures of g-priors for Bayesian Variable Selection. Journal of the
#' American Statistical Association. 103:410-423.  \cr
#' \url{https://dx.doi.org/10.1198/016214507000001337}
#' @keywords regression
#' @examples
#'
#' library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' # EB local uses a different g within each model
#' crime.EBL =  bas.lm(y ~ ., data=UScrime, n.models=2^15,
#'                     prior="EB-local", initprobs= "eplogp")
#' # use a common (global) estimate of g
#' crime.EBG = EB.global(crime.EBL)
#'
#'
#' @rdname EB.global
#' @family coef priors
#' @export
EB.global = function(object, tol= .1, g.0=NULL, max.iterations=100) {
n  = object$n
SSY = var(object$Y)*(n-1)
SSE <- (1.0 - object$R2)*SSY
SSR = SSY - SSE
p = object$size - 1
R2 = object$R2
prior = object$priorprobs

  postmodelprob <- function(R2, p, n, g, prior=1) {
     logmarg  <-  .5*((n - 1 - p)*log(1 + g) - (n-1)*log( 1 + g*(1 - R2)))
     logmarg[p == 0] = 0
     modelprob <- exp(logmarg -  max(logmarg))
     modelprob <- modelprob*prior/sum(modelprob*prior)
     if  (any(is.na(modelprob))) warning("NA's in modelprobs")
    return(modelprob)
   }

best = sort.list(-object$logmarg)[1]
sbest = min(object$shrinkage[best], .99)

if (is.null(g.0)) g.0 = sbest/(1 - sbest)
  tau.0 = g.0 + 1
  phi = (n - 1)/(SSY - (g.0/(1 + g.0))*SSR)
  post.prob = postmodelprob(object$R2,p, n, max(tau.0 - 1, 0), prior)
  tau.0 = sum(post.prob*SSR*phi)/(sum(post.prob*p))
  tau = tau.0 - 2*tol
  it = 0

while (abs(tau - tau.0) > tol | it < max.iterations) {
  g = max(tau - 1, 0)
  phi = (n - 1)/(SSY - (g/(1 + g))*SSR)
  post.prob = postmodelprob(object$R2,p, n, g, prior)
  tau.0 = tau
  tau = sum(post.prob*SSR*phi)/(sum(post.prob*p))
  it = it + 1
}

g = max(tau -1, 0)
logmarg  =  .5*((n -1 - p)*log(1 + g) - (n-1)*log( 1 + g*(1 - R2)))
logmarg[p == 0] = 0
postprobs = postmodelprob(object$R2,p, n, g, prior)
which = which.matrix(object$which, object$n.var)
object$probne0 = as.vector(postprobs %*% which)
object$postprobs=postprobs
object$g = g
object$logmarg = logmarg
object$shrinkage = object$shrinkage*0 + g/(1 + g)
object$method = "EB-global"
return(object)
}

