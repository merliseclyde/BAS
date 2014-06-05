\name{BAS-package}
\alias{BAS-package}
\alias{BAS}
\docType{package}
\title{Bayesian Model Averaging using Bayesian Adaptive Sampling}
\description{
Package for Bayesian Model Averaging in linear models
using stochastic or deterministic sampling without replacement
from posterior distributions. Prior distributions on
coefficients are of the form of  Zellner's g-prior or mixtures of g-priors.
Options include the Zellner-Siow Cauchy Priors, the Liang et
al hyper-g priors, Local and Global Empirical Bayes estimates of g, and
other default model selection criteria such as AIC
and BIC. Sampling probabilities may be updated based on the
sampled models.
}
\details{
\tabular{ll}{
Package: \tab BAS\cr
Version: \tab 0.9\cr
Date: \tab 2010-5-1\cr
Depends: \tab R (>= 2.8)\cr
License: \tab GPL-2\cr
URL: \tab http://www.stat.duke.edu/~clyde\cr
}

Index:
\preformatted{
}

}
\author{Merlise Clyde, 
  \cr
  Maintainer: Merlise Clyde <clyde@stat.duke.edu>
}
\references{
  Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling
  for Variable Selection and Model Averaging. Journal of Computational
  and Graphical Statistics. (to appear) Department of Statistical
  Science Discussion  Paper 2009-16. Duke University. \cr
  \url{}
  
  Clyde, M. and George, E. I. (2004) Model uncertainty. Statist. Sci.,
  19, 81-94. \cr \url{http://www.isds.duke.edu/~clyde/papers/statsci.pdf}

  Clyde, M. (1999)
  Bayesian Model Averaging and Model Search Strategies (with
  discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid,
  J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages
  157-185.

   Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger,
  J.O. (2008) Mixtures of  g-priors for Bayesian Variable
  Selection. Journal of the American Statistical Association. 103:410-423.  \cr \url{http://dx.doi.org/10.1198/016214507000001337}
}
\keyword{package}
\keyword{regression}
\seealso{
\code{\link[BAS]{bas}}
}
\examples{
demo(BAS.USCrime)
demo(BAS.hald)
}
