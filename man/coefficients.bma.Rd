\name{coef.bas}
\alias{coef.bas}
\alias{coef}
\alias{coefficients}
\alias{coefficients.bas}
\alias{print.coef.bas}
\title{Coefficients of a Bayesian Model Average object}

\description{Extract conditional posterior means and standard deviations,
  marginal posterior means and standard deviations, posterior
  probabilities, and marginal inclusions probabilities under
  Bayesian Model Averaging from an object of class 'bas'}
\usage{
\method{coef}{bas}(object, ...)
\method{print}{coef.bas}(x, n.models=5,digits = max(3, getOption("digits") - 3),...)
}
\arguments{
  \item{object}{object of class 'bas' created by BAS}
  \item{x}{object of class 'coef.bas' to print}
  \item{n.models}{Number of top models to report in the printed summary}
  \item{digits}{number of significant digits to print}
  \item{...}{other optional arguments }
}
\value{
  \code{coefficients} returns an object of class coef.bas with the following:
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
  Selection.  Journal of the American Statistical Association.
  103:410-423.  \cr \url{http://dx.doi.org/10.1198/016214507000001337}
 
}

\seealso{\code{\link{bas}}, \code{\link{confint.coef.bas}}
}
\examples{data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, n.models=2^4, alpha=13,
                      prior="ZS-null", initprobs="Uniform", update=10)
coef.hald.gprior = coefficients(hald.gprior)
coef.hald.gprior
plot(coef.hald.gprior)
}
\author{Merlise Clyde \email{clyde@stat.duke.edu} }
\note{With highly correlated variables,
  marginal summaries may not be representative of the joint
  distribution. Use \code{\link{plot.coef.bas}} to view distributions.}
\keyword{regression}
