\name{intrinsic}
\alias{intrinsic}

\title{Intrinsic Prior Distribution for Coefficients in BMA Models}
\description{Creates an object representing the intrinsic prior on g,
  a special case of the tCCH mixture of g-priors on
  coefficients for BAS. 
}

\usage(intrinsic(n)}


\arguments{\item{n}{the sample size}
 }

 \value{returns an object of class "prior", with the family "intrinsic"
    of class "TCCH"  and hyerparameters alpha = 1, 
        beta = 1, s = 0, r = 1, n = n) for the tCCH prior where
    theta in the tCCH prior is determined by the model size and sample size.
}
\details{Creates a structure used for \code{\link{bas.glm}}. 
    }
\seealso{\code{\link{tCCH}}, \code{\link{robust}}, \code{\link{hyper.g}},
  \code{\link{hyper.g.n}}\code{\link{bas.glm}}}
\examples{n = 500;
 tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n)
}

\references{
  Womack, A., Novelo,L.L., Casella, G. (2014). Inference From Intrinsic
  Bayes’ Procedures Under Model Selection and Uncertainty. JASA
  109:1040-1053.
  \url{http://amstat.tandfonline.com/doi/abs/10.1080/01621459.2014.880348} 

  }
\author{Merlise Clyde}

