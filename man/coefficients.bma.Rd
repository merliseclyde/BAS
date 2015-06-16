\name{coef.bma}
\alias{coef.bma}
\alias{coef}
\alias{coefficients}
\alias{coefficients.bma}
\alias{print.coef.bma}
\title{Coefficients of a Bayesian Model Average object}

\description{Extract conditional posterior means and standard deviations,
  marginal posterior means and standard deviations, posterior
  probabilities, and marginal inclusions probabilities under
  Bayesian Model Averaging from an object of class BMA}
\usage{
\method{coef}{bma}(object, ...)
\method{print}{coef.bma}(x, n.models=5,digits = max(3, getOption("digits") - 3),...)
}
\arguments{
  \item{object}{object of class 'bma' created by BAS}
  \item{x}{object of class 'coef.bma' to print}
  \item{n.models}{Number of top models to report in the printed summary}
  \item{digits}{number of significant digits to print}
  \item{...}{other optional arguments }
}
\value{
  \code{coefficients} returns an object of class coef.bma with the following:
  \item{conditionalmeans}{a  matrix with conditional posterior means
  for each model}
\item{conditionalsd}{ standard deviations for each model }
\item{postmean}{marginal posterior means of each regression coefficient
  using BMA}
\item{postsd}{marginal posterior standard deviations using BMA}
\item{postne0}{vector of posterior inclusion probabilities, marginal
  probability that a coefficient is non-zero}
}

\details{Calculates posterior means and (approximate) standard
  deviations of the regression coefficients under Bayesian Model
  averaging using g-priors and mixtures of g-priors.  Print returns
  overall summaries. For fully Bayesian
  methods that place a prior on g, the posterior standard deviations do
  not take into account full uncertainty regarding g. Will be updated in
  future releases.
 }

 \references{
 Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger,
  J.O. (2005) Mixtures of  g-priors for Bayesian Variable
  Selection.  \cr \url{http://www.stat.duke.edu/05-12.pdf}
  
}

\seealso{\code{\link{bas}}
}
\examples{data("Hald")
\dontrun{hald.gprior =  bas.lm(Y~ ., data=Hald, n.models=2^4, alpha=13,
                      prior="ZS-null", initprobs="Uniform", update=10)
coef.hald.gprior = coefficients(hald.gprior)
coef.hald.gprior
plot(coef.hald.gprior)
}}
\author{Merlise Clyde \email{clyde@stat.duke.edu} }
\note{With highly correlated variables,
  marginal summaries may not be representative of the
  distribution. Use \code{\link{plot.coef.bma}} to view distributions.}
\keyword{regression}
