\name{EB.global}
\alias{EB.global}
\alias{EB.global.bma}
\title{Finds the global Empirical Bayes estimates for BMA }
\description{Finds the  global Empirical Bayes estimates of g in
  Zellner's g-prior and model probabilities}
\usage{EB.global.bma(object, tol= .1, g.0=NULL, max.iterations=100)
}
\arguments{
  \item{object}{A 'bma' object created by \code{\link{bas}}}
  \item{tol}{tolerance for estimating g}
  \item{g.0}{intial value for g}
  \item{max.iterations}{Maximum number of iterations for the EM algorithm}
  }
\value{An object of class 'bma' using Zellner's g prior with an estimate
  of g based on all models}
\details{Uses the EM algorithm in Liang et al to estimate the type II
  MLE of g in Zellner's g prior} 
\references{  Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger,
  J.O. (2008) Mixtures of  g-priors for Bayesian Variable
  Selection. Journal of the American Statistical
  Association. 103:410-423.  \cr
  \url{http://dx.doi.org/10.1198/016214507000001337}
}
\seealso{\code{\link{bas}}, \code{\link{update}}}

\examples{
\dontrun{library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
# EB local uses a different g within each model
crime.EBL =  bas.lm(y ~ ., data=UScrime, n.models=2^15,
                    prior="EB-local", initprobs= "eplogp")
# use a common (global) estimate of g
crime.EBG = EB.global.bma(crime.EBL)
}}
\author{Merlise Clyde \email{clyde@stat.duke.edu}}
\keyword{regression}
