\name{bayesglm.fit}
\alias{bayesglm.fit}
\title{Fitting Generalized Linear Models Bayesian marginal likelihood evaluation}
\description{A version of glm.fit rewritten in C; also returns marginal
  likelihoods for Baysesian model comparison}
\usage{bayesglm.fit(x, y, weights = rep(1, nobs),
            start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, nobs), family = binomial(),
            coefprior = bic.prior(nobs),
            control = glm.control(), intercept = TRUE) 
}
\arguments{\item{x}{design matrix}
  \item{y}{response}
  \item{weights}{optional vector of weights to be used in the fitting
    process.  SHould be NULL or a numeric vector.}
  \item{start}{starting value for coefficients in the linear predictor}
  \item{etastart}{starting values for the linear predictor}
  \item{mustart}{starting values for the vectors of means}
  \item{offset}{a priori known component to be included in the linear
predictor}
  \item{family}{a description of the error distribution and link 
    function for exponential family;
    currently only binomial() is coded.}
  \item{coefprior}{function specifying prior distribution on
    coefficients with optionlal hyperparamters leading to marginal
    likelihood calculations; options  include \code{bic.prior()},\code{
      aic.prior()}, and     \code{ic.prior()}}
  \item{control}{a list of parameters that control convergence in the
    fitting process.  See the documentation for
    \code{glm.control()}}
  \item{intercept}{should an intercept be included in the null model?}
}
\value{
  \item{coefficients}{MLEs}
  \item{se}{Standard errors of coefficients based on the sqrt of the
    diagonal of the inverse information matrix}
  \item{mu}{fitted mean}
  \item{rank}{numeric rank of the fitted linear model}
  \item{deviance}{minus twice the log likelihood evaluated at
    the MLEs}
  \item{g}{value of g in g-priors}
  \item{shrinkage}{shrinkage factor for coefficients in linear
  predictor}
 \item{RegSS}{quadratic form beta'I(beta)beta used in shrinkage}
  \item{logmarglik}{the log marginal or integrated log likelihood (up to a constant)}
}
\details{C version of glm-fit.  For different prior choices returns,
  marginal likelihood  of model using a Laplace approximation.
}
\references{\code{\link{glm}}
}
\seealso{\code{\link{bic.prior}}
}
\examples{
require(MASS)
library(MASS)
Pima.tr
Y = as.numeric(Pima.tr$type) - 1
X = cbind(1, as.matrix(Pima.tr[,1:7]))
out = bayesglm.fit(X, Y, family=binomial(),coefprior=bic.prior(n=length(Y)))
out$coef
# using built in function
glm(type ~ ., family=binomial(), data=Pima.tr)
}
\author{Merlise Clyde translated the \code{\link{glm.fit}} from R base
        into C using the .Call interface}
\keyword{regression}
