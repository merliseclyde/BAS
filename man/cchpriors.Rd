\name{CCH}
\alias{CCH}

\title{Generalized g-Prior Distribution for Coefficients in BMA Models}
\description{Creates an object representing the CCH mixture of g-priors on
  coefficients for BAS . 
}

\usage{CCH(alpha, beta, s)}


\arguments{\item{alpha}{a scalar > 0, recommended alpha=.5 (betaprime)
    or 1 for CCH.  The hyper.g(alpha) is equivalent to CCH(alpha -2, 2, 0). Liang
    et al recommended values in the range 2 < alpha_h <= 4}
\item{beta}{a scalar > 0.  The value is not updated by the data; beta
  should be a function of n for consistency under the null model.  The
  hyper-g corresonds to b = 2}
\item{s}{a scalar, recommended s=0}
 
}

\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{IC.prior}}, \code{\link{bic.prior}}, \code{\link{bas.glm}}}
\examples{CCH(alpha=.5, beta=100, s=0) } 

\author{Merlise A Clyde}

