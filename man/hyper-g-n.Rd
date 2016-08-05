\name{hyper.g.n}
\alias{hyper.g.n}

\title{Generalized hyper-g/n Prior Distribution for g for mixtures of
  g-priors on  Coefficients in BMA Models}
\description{Creates an object representing the hyper-g/n mixture of g-priors on
  coefficients for BAS. This is a special case of the tCCH prior
}

\usage{hyper.g.n(alpha=3, n=NULL)}


\arguments{\item{alpha}{a scalar > 0, recommended 2 < alpha <= 3 }
\item{n}{The sample size; if NULL, the value derived from the data in the call to `bas.glm` will be used.}
}

\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}.  This is a
  special case of the \code{\link{tCCH}}, where \code{hyper.g.n(alpha=3,
  n)} is equivalent to \code{ tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n) } 
    }
\seealso{\code{\link{tCCH}}, \code{\link{robust}}, \code{\link{hyper.g}},
  \code{\link{CCH}}\code{\link{bas.glm}}}
\examples{n = 500
 hyper.g.n(alpha = 3, n=n)}
\author{Merlise Clyde}

