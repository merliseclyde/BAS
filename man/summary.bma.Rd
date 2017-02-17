\name{summary.bas}
\alias{summary.bas}
\alias{summary}
\alias{print.bas}
\alias{print}

\title{Summaries of Bayesian Model Averaging objects from BAS}
\description{
     \code{summary} and \code{print} methods for Bayesian model
     averaging objects created by \code{bas} Bayesian Adaptive Sampling
}

\usage{
\method{summary}{bas}(object, n.models = 5, ...)
\method{print}{bas}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{object of class 'bas'}
  \item{x}{object of class 'bas'}
  \item{n.models}{ optional number specifying the number of best models to display in summary }
  \item{digits}{ optional number specifying the number of digits to display }
%  \item{conditional}{ optional logical value specifying whether to display conditional expectation and standard deviation}
  \item{...}{other parameters to be passed to \code{print.default}}
}
\details{
The print methods display a view similar to \code{print.lm} . 
The summary methods display a view specific to Bayesian model averaging
giving the top 5 highest probability models represented by their inclusion indicators.  Summaries of the models include the Bayes Facotr (BF) of each model to the model with the largest marginal likelihood, the posterior probabilty of the models, R2, dim (which includes the intercept) and the log of the marginal likelihood.
}
\author{Merlise Clyde \email{clyde@stat.duke.edu}}
\seealso{\code{\link{coef.bas}}}
\examples{
library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
print(crime.bic)
summary(crime.bic)
}
\keyword{print}
\keyword{regression}
