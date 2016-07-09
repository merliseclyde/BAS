\name{update.bas}
\alias{update}
\alias{update.bas}
\title{Update BMA object using a new prior}
\description{Update a BMA object using a new prior distribution on the coefficients.
}
\usage{\method{update}{bas}(object, newprior, alpha=NULL, ...)
}
\arguments{
  \item{object}{BMA object to update}
  \item{newprior}{Update posterior model probabilities, probne0,
    shrinkage, logmarg, etc, using prior based on newprior.  See
    \code{\link{bas}} for available methods}
  \item{alpha}{optional new value of hyperparameter in prior for method}
  \item{...}{optional arguments}
}
\value{A new object of class BMA}
\details{Recomputes the marginal likelihoods for the new methods for
  models already sampled in current object. 
}
\references{
   Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling
  for Variable Selection and Model Averaging. Journal of Computational
  Graphics and Statistics.  20:80-101 \cr
  \url{http://dx.doi.org/10.1198/jcgs.2010.09049}
}
\seealso{\code{\link{bas}} for available methods and choices of alpha }
\examples{
\dontrun{
library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
crime.bic =  bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC",initprobs= "eplogp")
crime.zs = update(crime.bic, newprior="ZS-null")
}}
\author{Merlise Clyde \email{clyde@stat.duke.edu}}
\keyword{regression}
