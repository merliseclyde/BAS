\name{Jeffreys}
\alias{Jeffreys}

\title{Jeffreys Prior Distribution for $g$ for Mixtures of g-Priors for
  Coefficients in BMA Models} \description{Creates an object
  representing the Jeffrey's Prior on g mixture of g-priors on
  coefficients for BAS. This is equivalent to a limiting version of the
  CCH(a, 2, 0) with a = 0 or they hyper-g(a = 2) and is an improper
  prior.  As $g$ does not appear in the Null Model, Bayes Factors and
  model probabilities are not well-defined because of arbitrary
  normalizing constants, and for this reason the null model is excluded
  and the same c onstants are used across other models.  }
  \usage{Jeffreys()}


\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{CCH}} \code{\link{bas.glm}}}
\examples{Jeffreys() } 

\author{Merlise Clyde}

