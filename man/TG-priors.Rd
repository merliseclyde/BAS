\name{TG}
\alias{TG}

\title{Generalized g-Prior Distribution for Coefficients in BMA Models}
\description{Creates an object representing the Truncated Gamma (tCCH) mixture of g-priors on
  coefficients for BAS, where u = 1/(1+g) has a Gamma distribution
  supported on (0, 1]. 
}

\usage{TG(alpha)}


\arguments{\item{alpha}{a scalar > 0, recommended alpha=.5 (betaprime)
    or 1.  alpha=2 corresponds to the uniform prior on the shrinkage factor.}
}

\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{CCH}} \code{\link{bas.glm}}}
\examples{CCH(alpha=.5, beta=100, s=0) } 

\author{Merlise Clyde}

