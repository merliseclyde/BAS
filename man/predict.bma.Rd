\name{predict.bas}
\alias{predict.bas}
\alias{predict}
\title{Prediction Method for an object of class BMA}
\description{Predictions under model averaging or other estimators from
  a BMA object of class inheriting from 'bas'.
}
\usage{\method{predict}{bas}(object, newdata, se.fit=FALSE, type="link", top=NULL,
 estimator="BMA", prediction=FALSE, ...)}
\arguments{
  \item{object}{An object of class BAS, created by \code{bas}}
  \item{newdata}{dataframe, new matrix or vector of data for predictions. If missing, then use the data used for fitting for obtaining fitted and predicted values. May include a column for the intercept or just the predictor variables. If a dataframe, the variables are extracted using model.matrix using the call that created 'object'.  }
  \item{top}{a scalar integer M.  If supplied, subset the top M models,
  based on posterior probabilities for model predictions and BMA.} 
  \item{type}{Type of predictions required. "Link" which is on the
  scale of the linear predictor is the only option currently.}
 \item{estimator}{estimator used for predictions.  Currently supported
   options include:
    \cr
    'HPM' the highest probability model \cr
    'BMA' Bayesian model averaging, using optionally only the 'top'
    models \cr
    'MPM' the median probability model of Barbieri and Berger. \cr
    'BPM' the model that is closest to BMA predictions under squared
  error loss. BMA may be computed using only the 'top' models if supplied}
  \item{se.fit}{indicator for whether to compute se of fitted and predictied values}
  \item{prediction}{logical value to indicate whether the observed design matrix used in fitting or the newdata will be used for predictions}
  \item{...}{optional extra arguments}
}
\value{a list of
  \item{fit}{fitted values based on the selected estimator}
  \item{Ybma}{predictions using BMA, the same as fit for non-BMA methods for compatibilty}
  \item{Ypred}{matrix of predictions under each model for BMA}
  \item{se.fit}{se of fitted values}
  \item{se.pred}{se for predicted values}
  \item{best}{index of top models included}
  \item{bestmodels}{subset of bestmodels used for fitting or prediction}
}
\details{Use BMA and/or model selection to form predictions using the top highest probability models.}
\seealso{\code{\link{bas}}, \code{\link{fitted.bas}}
}
\examples{
data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
predict(hald.gprior, top=5, type="response", predict=FALSE, se.fit=TRUE)  # default is BMA
predict(hald.gprior, newdata=Hald, estimator="MPM", se.fit=TRUE, prediction=FALSE)
# same as fitted
fitted(hald.gprior,estimator="MPM")
}
\author{Merlise Clyde}
\keyword{regression}

