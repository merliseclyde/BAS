\name{testBF.prior}
\alias{testBF.prior}
\title{Test based Bayes Factors for BMA Models}
\description{Creates an object representing the prior distribution on
  coefficients for BAS that coreesponds to the test-based Bayes Factors. 
}
\usage{testBF.prior(g)
}
\arguments{\item{g}{a scalar used in the covariance of Zellner's
    g-prior, Cov(beta) = sigma^2 g (X'X)^-
    
    }
 }   
\value{returns an object of class "prior", with the family and
  hyerparameters.
}
\details{Creates a  prior object structure used for BAS in `bas.glm`. 
    }
\seealso{\code{\link{g.prior}}, \code{\link{bas.glm}}}
\examples{
testBF.prior(100)
library(MASS)
data(Pima.tr)

# use g = n 
bas.glm(type ~ ., data=Pima.tr, family=binomial(), 
        betaprior=testBF.prior(nrow(Pima.tr)),
        modelprior=uniform(), method="BAS")
} 

\author{Merlise Clyde}

