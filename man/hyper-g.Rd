\name{hyper.g}
\alias{hyper.g}

\title{Hyper-g-Prior Distribution for Coefficients in BMA Models}
\description{Creates an object representing the hyper-g mixture of g-priors on
  coefficients for BAS. 
}
\usage{hyper.g(alpha)}


\arguments{\item{alpha}{a scalar > 0. The hyper.g(alpha) is equivalent
    to CCH(alpha -2, 2, 0). Liang     et al recommended values in the
    range 2 < alpha_h <= 3}. 
}

\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{CCH}} \code{\link{bas.glm}}}
\examples{hyper.g(alpha=.5) } 

\author{Merlise Clyde}

