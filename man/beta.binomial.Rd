\name{beta.binomial}
\alias{beta.binomial}
\alias{Beta.Binomial}
\title{Beta-Binomial Prior Distribution for Models}
\description{Creates an object representing the prior distribution on
  models for BAS. 
}
\usage{beta.binomial(alpha=1.0, beta=1.0)
}
\arguments{\item{alpha}{parameter in the beta prior distribution}
\item{beta}{parameter in the beta prior distribution}
}
\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{The beta-binomial distribution on model size is obtained by
  assigning each variable inclusion indicator independent Bernoulli
  distributions with probability w, and then giving w a beta(alpha,beta)
  distribution. Marginalizing over w leads to the distribution on model
  size having the beta-binomial distribution. The default hyperparaeters
  lead to a uniform distribution over model size.   
}
\seealso{\code{\link{bas.lm}}, \code{\link{Bernoulli}},\code{\link{uniform}}
}
\examples{beta.binomial(1,10)} 

\author{Merlise Clyde}

