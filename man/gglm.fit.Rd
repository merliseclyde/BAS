\name{gglm.fit}
\alias{gglm.fit}
\title{
Fitting Generalized Linear Models.
}
\description{
Fitting Generalized Linear Models.A version of glm.fit rewritten in C.
}
\usage{
gglm.fit(x, y, weights = rep(1, nobs), offset = rep(0, nobs), 
 family = binomial(), control = glm.control())
}
\arguments{
  \item{x}{design matrix}
  \item{y}{response}
  \item{weights}{optional vector of weights to be used in the fitting
    process.  SHould be NULL or a numeric vector.}
  \item{offset}{a priori known component to be included in the linear
predictor}
  \item{family}{a description of the error distribution and link 
    function for exponential family;
    currently only binomial() is coded.}
  \item{control}{a list of parameters that control convergence in the
    fitting process.  See the documentation for
    \code{glm.control()}}
}
\value{
  \item{coefficients}{MLEs}
  \item{se}{Standard errors of coefficients based on the sqrt of the
    diagonal of the inverse information matrix}
  \item{mu}{fitted mean}
  \item{rank}{numeric rank of the fitted linear model}
  \item{deviance}{minus twice the log likelihood evaluated at
    the MLEs}
 \item{RegSS}{quadratic form beta'I(beta)beta used in shrinkage}
 \item{fitted.values}{fitted-values}
}
\examples{
library(BAS)
n <- 100
x <- 1 + rnorm(n,1) 
xbeta <- -1  + (x* 2)
proba <- exp(xbeta)/(1 + exp(xbeta))
y <- ifelse(runif(n,1) < proba,1,0)
#fit <- gglm.fit(cbind(1,x),y)
}
