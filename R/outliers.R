########## Functions for Chaloner and Brant ###############
#' Bayesian Outlier Detection
#'
#' Calculate the posterior  probability that the absolute value of
#' error exceeds more than k standard deviations
#' P(|epsilon_j| > k sigma | data)
#' under the model Y = X B + epsilon,
#' with epsilon ~ N(0, sigma^2 I)
#' based on the paper
#' by Chaloner & Brant Biometrika (1988). Either k or the prior
#' probability of there being no outliers must be provided.
#' This only uses the reference prior p(B, sigma) = 1;
#' other priors and model averaging to come.

#'
#' @param lmobj  An object of class `lm`
#' @param k number of standard deviations used in calculating
#' probability of an individual case being an outlier,
#' P(|error| > k sigma | data)
#' @param prior.prob The prior probability of there being no
#' outliers in the sample of size n
#' @return Returns a list of three items:
#' \item{e}{residuals}
#' \item{hat}{leverage values}
#' \item{prob.outlier}{posterior probabilities of a point being an outlier}
#' \item{prior.prob}{prior probability of a point being an outlier}
#' @references Chaloner & Brant (1988)
#' A Bayesian Approach to Outlier Detection and Residual Analysis
#' Biometrika (1988) 75, 651-659
#' @examples
#' data("stackloss")
#' stack.lm <- lm(stack.loss ~ ., data = stackloss)
#' stack.outliers <- Bayes.outlier(stack.lm, k = 3)
#' plot(stack.outliers$prob.outlier, type = "h", ylab = "Posterior Probability")
#' # adjust for sample size for calculating prior prob that a
#' # a case is an outlier
#' stack.outliers <- Bayes.outlier(stack.lm, prior.prob = 0.95)
#' # cases where posterior probability exceeds prior probability
#' which(stack.outliers$prob.outlier > stack.outliers$prior.prob)
#' @export
Bayes.outlier <- function(lmobj, k, prior.prob) {
  e <- residuals(lmobj)
  h <- hatvalues(lmobj)
  alpha <- (lmobj$df.residual) / 2
  rate <- (lmobj$df.residual * (summary(lmobj)$sigma)^2) / 2
  n <- length(e)

  if (missing(k) & missing(prior.prob)) {
    stop("please provide either k or the prior probability of no outliers")
  }
  else {
    if (missing(k)) k <- qnorm(.5 + .5 * (prior.prob^(1 / n)))
  }
  pr <- rep(0, n)
  for (i in 1:n) {
    pr[i] <- integrate(
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
    prob.outlier = pr,
    prior.prob = pnorm(-k) * 2
  ))
}

outlier.prob <- function(phi, ehat, hii, alpha, rate, nsd) {
  z1 <- (nsd - ehat * sqrt(phi)) / sqrt(hii)
  z2 <- (-nsd - ehat * sqrt(phi)) / sqrt(hii)
  pr.phi <- (1 - pnorm(z1) + pnorm(z2)) * dgamma(phi, shape = alpha, rate = rate)
  return(pr.phi)
}
