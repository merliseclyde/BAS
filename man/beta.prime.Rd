\name{beta.prime}
\alias{beta.prime}
\title{ Beta-Prime Prior Distribution for Coefficients in BMA Model}
\description{Creates an object representing the Beta-Prime prior  that is  mixture of g-priors on coefficients for BAS. 
}
\usage{beta.prime(n)
}
\arguments{\item{n}{the sample size.}
 }

\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{CCH}}}
\examples{beta.prime(n=100)} 

\author{Merlise Clyde}

