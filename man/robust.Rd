\name{robust}
\alias{robust}
\title{ Robust-Prior Distribution for Coefficients in BMA Model}
\description{Creates an object representing the robust prior of Bayarri
  et al (2012) that is  mixture of g-priors on coefficients for BAS. 
}
\usage{robust(n)
}
\arguments{\item{n}{the sample size.}
 }

\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{CCH}} and\code{\link{bas.glm}} }
\examples{robust(100)} 

\author{Merlise Clyde}

