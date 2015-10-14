\name{cv.summary.bma}
\alias{cv.summary.bma}

\title{Summaries for Out of Sample Prediction}

\description{Compute summaries from out of sample predictions for a BMA object}

\usage{cv.summary.bma(object, pred, ytrue)}

\arguments{
  \item{object}{am object of class 'bma'}
  \item{pred}{output from \code{\link{predict.bas}}}
  \item{ytrue}{vector of left out response values}
}

\value{A matrix with the best models, posterior probabilities, R2,
  dimensions, Average Prediction Error from the HPM and Average
  prediction error for BMA prediction}

 
\seealso{\code{\link{predict.bas}}}

\author{Merlise Clyde \email{clyde@stat.duke.edu}}
\keyword{regression}
