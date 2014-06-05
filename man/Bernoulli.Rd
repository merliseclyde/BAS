\name{Bernoulli}
\alias{bernoulli}
\alias{Bernoulli}
\title{Independent Bernoulli Prior Distribution for Models}
\description{Creates an object representing the prior distribution on
  models for BAS. 
}
\usage{Bernoulli(probs=.5)
}
\arguments{\item{probs}{a scalar or vector of prior inclusion
    probabilities. If a scalar, the values is replicated for all
    variables ans a 1 is added for the intercept. BAS checs to see if
    the length is equal to the dimension of the parameter vector for the
  full model and addas a 1 to include the intercept.}
 }   
\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{The independent Bernoulli prior distribution is a commonly used
  prior in BMA, with the Uniform distribution a special case with probs=.5.  If all
  indicator variables have a independent Bernoulli
  distributions with common probability probs, the distribution on model
  size binomial(p, probs) distribution.
  }
\seealso{\code{\link{bas.lm}}, \code{\link{beta.binomial}},\code{\link{uniform}
}}
\examples{Bernoulli(.9)} 

\author{Merlise Clyde}

