\name{summary.bma}
\alias{summary.bma}
\alias{summary}
\alias{print.bma}
\alias{print}

\title{Summaries of Bayesian Model Averaging objects}
\description{
     \code{summary} and \code{print} methods for Bayesian model
     averaging objects created by \code{bas} Bayesian Adaptive Sampling
}

\usage{
\method{summary}{bma}(object, n.models = 5, ...)
\method{print}{bma}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{object of class 'bma'}
  \item{x}{object of class 'bma'}
  \item{n.models}{ optional number specifying the number of best models to display in summary }
  \item{digits}{ optional number specifying the number of digits to display }
%  \item{conditional}{ optional logical value specifying whether to display conditional expectation and standard deviation}
  \item{...}{other parameters to be passed to \code{print.default}}
}
\details{
The print methods display a view similar to \code{print.lm} . 
The summary methods display a view specific to Bayesian model averaging
giving the top highest probability models.
}
\author{Merlise Clyde \email{clyde@stat.duke.edu}}
\seealso{\code{\link{coefficients.bma}}}
\examples{
\dontrun{library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
print(crime.bic)
summary(crime.bic)
}}
\keyword{print}
\keyword{regression}
