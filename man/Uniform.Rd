\name{uniform}
\alias{uniform}
\alias{Uniform}
\title{Uniform Prior Distribution for Models}
\description{Creates an object representing the prior distribution on
  models for BAS. 
}
\usage{uniform()
}

\value{returns an object of class "prior", with the family name Uniform.
  }
\details{The Uniform  prior distribution is a commonly used
  prior in BMA, and is a special case of the indpendent Bernoulli prior
  with probs=.5. The implied prior distribution on model
  size is binomial(p, .5).
  }
\seealso{\code{\link{bas.lm}}, \code{\link{beta.binomial}},\code{\link{Bernoulli}},
}
\examples{uniform()} 

\author{Merlise Clyde}

