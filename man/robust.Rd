% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_priors.R
\name{robust}
\alias{robust}
\title{Robust-Prior Distribution for Coefficients in BMA Model}
\usage{
robust(n = NULL)
}
\arguments{
\item{n}{the sample size.}
}
\value{
returns an object of class "prior", with the family and
hyerparameters.
}
\description{
Creates an object representing the robust prior of Bayarri et al (2012) that
is mixture of g-priors on coefficients for BAS.
}
\details{
Creates a prior structure used for \code{\link{bas.glm}}.
}
\examples{
robust(100)
}
\seealso{
\code{\link{CCH}} and\code{\link{bas.glm}}

Other beta priors: 
\code{\link{CCH}()},
\code{\link{EB.local}()},
\code{\link{IC.prior}()},
\code{\link{Jeffreys}()},
\code{\link{TG}()},
\code{\link{beta.prime}()},
\code{\link{g.prior}()},
\code{\link{hyper.g}()},
\code{\link{hyper.g.n}()},
\code{\link{intrinsic}()},
\code{\link{tCCH}()},
\code{\link{testBF.prior}()}
}
\author{
Merlise Clyde
}
\concept{beta priors}
