\name{tCCH}
\alias{tCCH}

\title{Generalized tCCH g-Prior Distribution for Coefficients in BMA Models}
\description{Creates an object representing the tCCH mixture of g-priors on
  coefficients for BAS. 
}

\usage{tCCH(alpha, beta, s, r, v, theta)}


\arguments{\item{alpha}{a scalar > 0, recommended alpha=.5 (betaprime)
    or 1.}
\item{beta}{a scalar > 0.  The value is not updated by the data; beta
  should be a function of n for consistency under the null model.}
\item{s}{a scalar, recommended s=0 a priori}
\item{r}{r arbitrary; in the hyper-g-n prior sets r = (alpha + 2)}
\item{v}{0 < v }{ 1/v determines the upper bound on the support of 1/(1
    + g). }
\item{theta}{theta > 1} {recommend theta = 1/n as inn the hyper-g-n}
 }

\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{CCH}}, \code{\link{robust}}, \code{\link{hyper.g}},
  \code{\link{hyper.g.n}}\code{\link{bas.glm}}}
\examples{n = 500;
 tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n)
}

\author{Merlise Clyde}

