\name{update.bma}
\alias{update}
\alias{update.bma}
\title{Update BMA object using a new prior}
\description{Update a BMA object using a new prior distribution on the coefficients.
}
\usage{\method{update}{bma}(object, newprior, alpha=NULL, ...)
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
  Clyde, M. and George, E. I. (2004) Model uncertainty. Statist. Sci.,
  19, 81-94. \cr
  \url{http://www.isds.duke.edu/~clyde/papers/statsci.pdf}
  \cr
  
  Clyde, M. (1999)
  Bayesian Model Averaging and Model Search Strategies (with
  discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid,
  J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages
  157-185.
  \cr
  
  Hoeting, J. A., Madigan, D., Raftery, A. E. and Volinsky, C. T. (1999)
  Bayesian model averaging: a tutorial (with discussion). Statist. Sci.,
  14, 382-401. \cr
  \url{http://www.stat.washington.edu/www/research/online/hoeting1999.pdf}
  \cr
  
  Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger,
  J.O. (2008) Mixtures of  g-priors for Bayesian Variable
  Selection.  JASA \cr
  \url{http://www.stat.duke.edu/05-12.pdf}
  \cr
  
  Zellner, A. (1986) On assessing prior distributions and Bayesian
  regression analysis with g-prior distributions. In Bayesian Inference
  and Decision Techniques: Essays in Honor of Bruno de Finetti,
  pp. 233-243. North-Holland/Elsevier.
  \cr
  
  Zellner, A. and Siow, A. (1980) Posterior odds ratios for selected
  regression hypotheses. In Bayesian Statistics: Proceedings of the First
  International Meeting held in Valencia (Spain), pp. 585-603.
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
