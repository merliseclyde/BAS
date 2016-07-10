\name{confint.coef.bas}
\alias{confint.coef.bas}
\alias{confint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Credible Intervals for BMA coefficients}
\description{Computes Credible Intervals for regression coefficients from BAS objects}

\usage{
\method{confint}{coef.bas}(object, parm, level = 0.95, nsim=10000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a coef.bas object
}
  \item{parm}{a specification of which parameters are to be given credible intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
  }
  \item{level}{the probabilty coverage required}
 
  \item{nsim}{number of Monte Carlo draws from the posterior distribution.  Used when number of models is greater than 1.}
  \item{...}{other arguments to passed; none currently}
}
\details{
Uses Monte Carlo simulations using posterior means and standard deviations of coefficents to generate draws from the posterior distributions and returns highest posterior density (HPD) credible intervals.  If the number of models equals one, then use the t distribution to find intervals.  These currently condition on the estimate of $g$.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix (or vector) with columns giving lower and upper HPD credible limits for each parameter. These will be labelled as  1-level)/2 and 1 - (1-level)/2 in percent
(by default 2.5 and 97.5).}



\author{Merlise A Clyde}

\note{For mixture of g-priors these are approximate.  This uses Monte Carlo sampling so results may be subjet to Monte Carlo variation and larger values of nsim may be needed to reduce variability.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
coef.hald = coef(hald.gprior)
confint(coef.hald)
confint(coef.hald, approx=FALSE, nsim=5000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")

