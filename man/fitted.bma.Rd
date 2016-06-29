\name{fitted.bas}
\alias{fitted.bas}
\alias{fitted}
\title{Fitted values for a BAS BMA objects}
\description{Calculate fitted values for a BAS BMA object}

\usage{\method{fitted}{bas}(object,  type="HPM", top=NULL, ...)}
\arguments{
  \item{object}{An object of class 'bma' as created by \code{\link{bas}}}
  \item{type}{type of fitted value to return.  Options include
   \cr
    'HPM' the highest probability model \cr
    'BMA' Bayesian model averaging, using optionally only the 'top'
    models \cr
    'MPM' the median probability model of Barbieri and Berger.
    'BPM' the model that is closest to BMA predictions under squared error loss}
\item{top}{optional argument specifying that the 'top' models will be
    used in constructing the BMA prediction, if NULL all models will be
    used.  If top=1, then this is equivalent to 'HPM'}
  \item{...}{optional arguments, not used currently}
}
\value{A vector of length n of fitted values.}
\details{Calcuates fitted values at observed design matrix using either
  the highest probability model, 'HPM', the posterior mean (under BMA)
  'BMA', the median probability model 'MPM' or the best predictive model 'BPM". 
  The median probability
  model is defined by including variable where the marginal inclusion
  probability is greater than or equal to 1/2. For type="BMA", the
  weighted average may be based on using a subset of the highest
  probability models if an optional argument is given for top.  By
  default BMA uses all sampled models, which may take a while to compute
  if the number of variables or number of models is large. 
  The "BPM" is found be computing the squared distance of the vector of fitted 
  values for a model and the fitted values under BMA and returns the model with the
  smallest  distance.  In the presence of multicollinearity this may be quite different from the MPM, with extreme collinearity may drop relevant predictors.
}
\references{
  Barbieri, M.  and Berger, J.O. (2004) Optimal predictive
  model selection. Annals of Statistics. 32, 870-897. \cr
  \url{http://projecteuclid.org/Dienst/UI/1.0/Summarize/euclid.aos/1085408489}
  
  Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling
  for Variable Selection and Model Averaging. Journal of Computational
  Graphics and Statistics.  20:80-101 \cr \url{http://dx.doi.org/10.1198/jcgs.2010.09049}  
}
\seealso{\code{\link{predict.bas}}
}
\examples{
data(Hald)
hald.gprior =  bas.lm(Y~ ., data=Hald, prior="ZS-null", initprobs="Uniform")
plot(Hald$Y, fitted(hald.gprior, type="HPM"))
plot(Hald$Y, fitted(hald.gprior, type="BMA"))
plot(Hald$Y, fitted(hald.gprior, type="MPM"))
plot(Hald$Y, fitted(hald.gprior, type="BPM"))
}
\author{Merlise Clyde \email{clyde@AT@stat.duke.edu}}
\keyword{regression}
