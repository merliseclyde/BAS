\name{cv.summary.bas}
\alias{cv.summary.bas}

\title{Summaries for Out of Sample Prediction}

\description{Compute average prediction error from out of sample predictions}

\usage{cv.summary.bas(pred, ytrue, score="squared-error")}

\arguments{
  \item{pred}{fitted or predicted value from the output from \code{\link{predict.bas}}}
  \item{ytrue}{vector of left out response values}
  \item{score}{function used to summarize error rate.  Either "squared-error", "percent-explained", or "miss-class"}
}

\value{For squared error, the average  prediction error for the Bayesian estimator
       error = sqrt(sum(ytrue - yhat)^2/npred) 
       while for binary data the missclassification rate is more appropriate.
       For continuous data the "percent-explained"  reports 
       ar, similar to an out of sample   R2. 
       }

 
\seealso{\code{\link{predict.bas}}}

\author{Merlise Clyde \email{clyde@stat.duke.edu}}

\examples{

\dontrun{
library(foreign)
cognitive = read.dta("http://www.stat.columbia.edu/~gelman/arm/examples/child.iq/kidiq.dta")
cognitive$mom_work = as.numeric(cognitive$mom_work > 1)
cognitive$mom_hs =  as.numeric(cognitive$mom_hs > 0)
colnames(cognitive) = c("kid_score", "hs","iq", "work", "age")

set.seed(42)
n = nrow(cognitive)
test = sample(1:n, size=round(.20*n), replace=FALSE)
testdata =  cognitive[test,]
traindata = cognitive[-test,]
cog_train = bas.lm(kid_score ~ ., prior="BIC", modelprior=uniform(), data=traindata)
yhat = predict(cog_train, newdata=testdata, estimator="BMA", se=F)
cv.summary.bas(yhat$fit, testdata$kid_score)
}

}
\keyword{regression}
