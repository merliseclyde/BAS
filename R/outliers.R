########## Functions for Chaloner and Brant ###############
#' Bayesian Outlier Detection
#'
#' Calculate the posterior  probability that the absolute value of
#' error exceeds more than k standard deviations
#' P(|epsilon_j| > k sigma | data)
#' under the model Y = X B + epsilon,
#' with epsilon ~ N(0, sigma^2 I)
#' based on the paper
#' by Chaloner & Brant Biometrika (1988).
#' This only uses the reference prior p(B, sigma) = 1;
#' other priors and model averaging to come.

#'
#' @param lmobj  An object of class `lm`
#' @param k number of standard deviations used in calcuating
#' default is `k=3`,
#' @return Returns a list of three items:
#' \item{e}{residuals}
#' \item{hat}{leverage values}
#' \item{pr}{posterior probabilities of a point being an outlier}
#' @references Chaloner & Brant (1988)
#' A Bayesian Approach to Outlier Detection and Residual Analysis
#' Biometrika (1988) 75, 651-659
#' @examples
#' library(MASS)
#' data(stackloss)
#' stack.lm =lm(stack.loss ~ ., data=stackloss)
#' stack.outliers = Bayes.outlier(stack.lm)
#' plot(stack.outliers$pr, type="h", ylab="Posterior Probability")
#' @export
Bayes.outlier <- function(lmobj, k = 3) {
  e <- residuals(lmobj)
  h <- hatvalues(lmobj)
  alpha <- (lmobj$df.residual) / 2
  rate <- (lmobj$df.residual * (summary(lmobj)$sigma) ^ 2) / 2
  n <- length(e)
  pr <- rep(0, n)
  for (i in 1:n) {
    pr[i] = integrate(
      outlier.prob,
      lower = 0,
      upper = Inf,
      ehat = e[i],
      hii = h[i],
      alpha = alpha,
      rate = rate,
      nsd = k
    )$value
  }
  return(list(
    e = e,
    hat = h,
    prob.outlier = pr
  ))
}

outlier.prob <- function(phi, ehat,hii,alpha,rate, nsd) {
	z1 <- (nsd - ehat*sqrt(phi))/sqrt(hii)
	z2 <- (- nsd - ehat*sqrt(phi))/sqrt(hii)
	pr.phi <- (1 - pnorm(z1) + pnorm(z2))*dgamma(phi,shape=alpha, rate=rate)
	return(pr.phi)}



