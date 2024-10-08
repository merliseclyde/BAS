% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_priors.R
\name{g.prior}
\alias{g.prior}
\title{Families of G-Prior Distribution for Coefficients in BMA Models}
\usage{
g.prior(g)
}
\arguments{
\item{g}{a scalar used in the covariance of Zellner's g-prior, Cov(beta) =
sigma^2 g (X'X)^-1}
}
\value{
returns an object of class "prior", with the family and
hyerparameters.
}
\description{
Creates an object representing the g-prior distribution on coefficients for
BAS.
}
\details{
Creates a structure used for BAS.
}
\examples{
g.prior(100)
}
\seealso{
\code{\link{IC.prior}}

Other beta priors: 
\code{\link{CCH}()},
\code{\link{EB.local}()},
\code{\link{IC.prior}()},
\code{\link{Jeffreys}()},
\code{\link{TG}()},
\code{\link{beta.prime}()},
\code{\link{hyper.g}()},
\code{\link{hyper.g.n}()},
\code{\link{intrinsic}()},
\code{\link{robust}()},
\code{\link{tCCH}()},
\code{\link{testBF.prior}()}
}
\author{
Merlise Clyde
}
\concept{beta priors}
