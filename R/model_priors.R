#' Uniform Prior Distribution for Models
#'
#' Creates an object representing the prior distribution on models for BAS.
#'
#' The Uniform prior distribution is a commonly used prior in BMA, and is a
#' special case of the independent Bernoulli prior with probs=.5. The implied
#' prior distribution on model size is binomial(p, .5).
#'
#' @aliases uniform Uniform
#' @return returns an object of class "prior", with the family name Uniform.
#' @author Merlise Clyde
#' @seealso \code{\link{bas.lm}},
#' \code{\link{beta.binomial}},\code{\link{Bernoulli}},
#' @examples
#' uniform()
#' @rdname uniform
#' @family priors modelpriors
#' @export
uniform <- function() {
  structure(list(family = "Uniform", hyper.parameters = .5), class = "prior")
}



#' Independent Bernoulli Prior Distribution for Models
#'
#' Creates an object representing the prior distribution on models for BAS.
#'
#' The independent Bernoulli prior distribution is a commonly used prior in
#' BMA, with the Uniform distribution a special case with probs=.5.  If all
#' indicator variables have a independent Bernoulli distributions with common
#' probability probs, the distribution on model size binomial(p, probs)
#' distribution.
#'
#' @aliases bernoulli Bernoulli
#' @param probs a scalar or vector of prior inclusion probabilities. If a
#' scalar, the values is replicated for all variables ans a 1 is added for the
#' intercept. BAS checks to see if the length is equal to the dimension of the
#' parameter vector for the full model and addas a 1 to include the intercept.
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{bas.lm}},
#' \code{\link{beta.binomial}},\code{\link{uniform} }
#' @examples
#' Bernoulli(.9)
#' @family priors modelpriors
#' @export
#' @rdname Bernoulli
Bernoulli <- function(probs = 0.5) {
  if (length(probs) == 1) {
    if (probs == .5) {
      structure(list(family = "Uniform", hyper.parameters = .5), class = "prior")
    } else {
      structure(list(family = "Bernoulli", hyper.parameters = probs), class = "prior")
    }
  }
  else {
    structure(list(family = "Bernoulli", hyper.parameters = probs), class = "prior")
  }
}



#' Beta-Binomial Prior Distribution for Models
#'
#' Creates an object representing the prior distribution on models for BAS.
#'
#' The beta-binomial distribution on model size is obtained by assigning each
#' variable inclusion indicator independent Bernoulli distributions with
#' probability w, and then giving w a beta(alpha,beta) distribution.
#' Marginalizing over w leads to the distribution on model size having the
#' beta-binomial distribution. The default hyperparaeters lead to a uniform
#' distribution over model size.
#'
#' @aliases beta.binomial Beta.Binomial
#' @param alpha parameter in the beta prior distribution
#' @param beta parameter in the beta prior distribution
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{bas.lm}}, \code{\link{Bernoulli}},\code{\link{uniform}}
#' @examples
#' beta.binomial(1, 10) #' @family priors modelpriors
#' @family priors modelpriors
#' @rdname beta.binomial
#' @export
beta.binomial <- function(alpha = 1.0, beta = 1.0) {
  structure(list(family = "Beta-Binomial", hyper.parameters = c(alpha, beta)),
    class = "prior"
  )
}




#' Truncated Beta-Binomial Prior Distribution for Models
#'
#' Creates an object representing the prior distribution on models for BAS
#' using a truncated Beta-Binomial Distribution on the Model Size
#'
#' The beta-binomial distribution on model size is obtained by assigning each
#' variable inclusion indicator independent Bernoulli distributions with
#' probability w, and then giving w a beta(alpha,beta) distribution.
#' Marginalizing over w leads to the number of included
#' predictors having a beta-binomial distribution. The default hyperparameters
#' lead to a uniform distribution over model size.  The Truncated version
#' assigns zero probability to all models of size > trunc.
#'
#' @aliases tr.beta.binomial tr.Beta.Binomial
#' @param alpha parameter in the beta prior distribution
#' @param beta parameter in the beta prior distribution
#' @param trunc parameter that determines truncation in the distribution i.e.
#' P(M; alpha, beta, trunc) = 0 if M > trunc.
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{bas.lm}}, \code{\link{Bernoulli}},\code{\link{uniform}}
#' @examples
#'
#' tr.beta.binomial(1, 10, 5)
#' library(MASS)
#' data(UScrime)
#' UScrime[, -2] <- log(UScrime[, -2])
#' crime.bic <- bas.lm(y ~ .,
#'   data = UScrime, n.models = 2^15, prior = "BIC",
#'   modelprior = tr.beta.binomial(1, 1, 8),
#'   initprobs = "eplogp"
#' )
#' @family priors modelpriors
#' @rdname tr.beta.binomial
#' @export
tr.beta.binomial <- function(alpha = 1.0, beta = 1.0, trunc) {
  structure(list(family = "Trunc-Beta-Binomial", hyper.parameters = c(alpha, beta, trunc)),
    class = "prior"
  )
}



#' Truncated Power Prior Distribution for Models
#'
#' Creates an object representing the prior distribution on models for BAS
#' using a truncated Distribution on the Model Size where the probability of
#' gamma = p^-kappa |gamma| where gamma is the vector of model indicators
#'
#' The beta-binomial distribution on model size is obtained by assigning each
#' variable inclusion indicator independent Bernoulli distributions with
#' probability w, and then giving w a beta(alpha,beta) distribution.
#' Marginalizing over w leads to the number of included
#' predictors having a beta-binomial distribution. The default hyperparameters
#' lead to a uniform distribution over model size.  The Truncated version
#' assigns zero probability to all models of size > trunc.
#'
#' @aliases tr.power.prior tr.Power.Prior
#' @param kappa parameter in the prior distribution that controls sparsity
#' @param trunc parameter that determines truncation in the distribution i.e.
#' P(gamma; alpha, beta, trunc) = 0 if |gamma| > trunc.
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{bas.lm}}, \code{\link{Bernoulli}},\code{\link{uniform}}
#' @examples
#'
#' tr.power.prior(2, 8)
#' library(MASS)
#' data(UScrime)
#' UScrime[, -2] <- log(UScrime[, -2])
#' crime.bic <- bas.lm(y ~ .,
#'   data = UScrime, n.models = 2^15, prior = "BIC",
#'   modelprior = tr.power.prior(2, 8),
#'   initprobs = "eplogp"
#' )
#' @family priors modelpriors
#' @rdname tr.power.prior
#' @export
tr.power.prior <- function(kappa = 2, trunc) {
  structure(list(family = "Trunc-Power-Prior", hyper.parameters = c(kappa, trunc)),
    class = "prior"
  )
}



#' Truncated Poisson Prior Distribution for Models
#'
#' Creates an object representing the prior distribution on models for BAS
#' using a truncated Poisson Distribution on the Model Size
#'
#' The Poisson prior distribution on model size is obtained by assigning each
#' variable inclusion indicator independent Bernoulli distributions with
#' probability w, and then taking a limit as p goes to infinity and w goes to
#' zero, such that p*w converges to lambda.  The Truncated version assigns zero
#' probability to all models of size M > trunc.
#'
#' @aliases tr.Poisson tr.poisson
#' @param lambda parameter in the Poisson distribution representing expected
#' model size with infinite predictors
#' @param trunc parameter that determines truncation in the distribution i.e.
#' P(M; lambda, trunc) = 0 if M > trunc
#' @return returns an object of class "prior", with the family and
#' hyerparameters.
#' @author Merlise Clyde
#' @seealso \code{\link{bas.lm}}, \code{\link{Bernoulli}},\code{\link{uniform}}
#' @examples
#' tr.poisson(10, 50)
#' @family priors modelpriors
#' @rdname tr.poisson
#' @export
tr.poisson <- function(lambda, trunc) {
  structure(list(family = "Trunc-Poisson", hyper.parameters = c(lambda, trunc)),
    class = "prior"
  )
}

#' Independent Bernoulli prior on models that with constraints for
#' model hierarchy induced by interactions
#' @param pi  Bernoulli probability that term is included
#' @param parents matrix of terms and parents with indicators of which terms
#' are parents for each term
#' @family priors modelpriors
#' @rdname Bernoulli.heredity
#' @note Not implemented yet for use with bas.lm or bas.glm
#' @export
Bernoulli.heredity <- function(pi = 0.5, parents) {
  stop("not implemented full yet")
  structure(list(
    family = "Bernoulli.Constrained",
    hyper.parameters = c(hyper.parameters = pi, parents = parents)
  ),
  class = "prior"
  )
}
