\name{IC.prior}
\alias{IC.prior}
\alias{aic.prior}
\alias{AIC.prior}
\alias{bic.prior}
\alias{BIC.prior}
\title{Information Criterion Families of Prior Distribution for Coefficients in BMA Models}
\description{Creates an object representing the prior distribution on
  coefficients for BAS. 
}
\usage{IC.prior(penalty)
}
\arguments{\item{penalty}{a scalar used in the penalized loglikelihood
    of the form penalty*dimension
    }
 }   
\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{The log marginal likelihood is approximated as -2*(deviance +
  penalty*dimension).  Allows alternatives to AIC (penalty = 2) and BIC
  (penalty = log(n)).  For BIC, the argument may be missing, in which case the sample size is determined from the call to `bas.glm` and used to dertermine the penalty.
    }
\seealso{\code{\link{g.prior}}}
\examples{IC.prior(2)
          aic.prior()
          bic.prior(100)} 

\author{Merlise Clyde}

