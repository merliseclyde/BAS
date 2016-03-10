\name{tr.beta.binomial}
\alias{tr.beta.binomial}
\alias{tr.Beta.Binomial}
\title{Truncated Beta-Binomial Prior Distribution for Models}
\description{Creates an object representing the prior distribution on
  models for BAS using a truncated Beta-Binomial Distribution on the
  Model Size 
}
\usage{tr.beta.binomial(alpha=1.0, beta=1.0, trunc)
}
\arguments{\item{alpha}{parameter in the beta prior distribution}
  \item{beta}{parameter in the beta prior distribution}
  \item{trunc}{parameter that determines truncation in the distribution
    i.e. P(M; alpha, beta, trunc) = 0 if  M > trunc}
}
\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{The beta-binomial distribution on model size is obtained by
  assigning each variable inclusion indicator independent Bernoulli
  distributions with probability w, and then giving w a beta(alpha,beta)
  distribution. Marginalizing over w leads to the distribution on model
  size having the beta-binomial distribution. The default hyperparaeters
  lead to a uniform distribution over model size.   The Truncated
  version assigns zero probability to all models of size > trunc.
}
\seealso{\code{\link{bas.lm}}, \code{\link{Bernoulli}},\code{\link{uniform}}
}
\examples{tr.beta.binomial(1,10, 25)} 

\author{Merlise Clyde}

