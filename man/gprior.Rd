\name{g.prior}
\alias{g.prior}
\title{Families of G-Prior Distribution for Coefficients in BMA Models}
\description{Creates an object representing the g-prior distribution on
  coefficients for BAS. 
}
\usage{g.prior(g)
}
\arguments{\item{g}{a scalar used in the covariance of Zellner's
    g-prior, Cov(beta) = sigma^2 g (X'X)^-1
    
    }
 }   
\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for BAS. 
    }
\seealso{\code{\link{IC.prior}}}
\examples{g.prior(100)} 

\author{Merlise Clyde}

