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
  \item{newdata}{dataframe for predictions. If missing, then use the dataframe used for fitting for obtaining fitted and predicted values.  }
   \item{se.fit}{indicator for whether to compute se of fitted and predictied values}
  \item{type}{Type of predictions required. "Link" which is on the
  scale of the linear predictor is the only option currently.}
  \item{top}{a scalar integer M.  If supplied, subset the top M models,
  based on posterior probabilities for model predictions and BMA.} 
 \item{estimator}{estimator used for predictions.  Currently supported
   options include:
    \cr
    'HPM' the highest probability model \cr
    'BMA' Bayesian model averaging, using optionally only the 'top'
    models \cr
    'MPM' the median probability model of Barbieri and Berger. \cr
    'BPM' the model that is closest to BMA predictions under squared
  error loss. BMA may be computed using only the 'top' models if supplied}
   \item{prediction}{logical value to indicate whether the observed design matrix used in fitting or the newdata will be used for estimating the mean or for making predictions.  The default is FALSE for estimation of the mean.}
  \item{...}{optional extra arguments}
}
\value{a list of
  \item{fit}{fitted values based on the selected estimator}
  \item{Ybma}{predictions using BMA, the same as fit for non-BMA methods for compatibilty; will be deprecated}
  \item{Ypred}{matrix of predictions under each model for BMA}
  \item{se.fit}{se of fitted values; in the case of BMA this will be a matrix}
  \item{se.pred}{se for predicted values; in the case of BMA this will be a matrix}
  \item{se.bma.fit}{vector of posterior sd under BMA for posterior mean of the 
  regression function. this will be NULL if estimator is not 'BMA'}
  \item{se.bma.pred}{vector of posterior sd under BMA for posterior predictive values.
    this will be NULL if estimator is not 'BMA'}
  \item{best}{index of top models included}
  \item{bestmodels}{subset of bestmodels used for fitting or prediction}
  \item{df}{scalar or vector of degrees of freedom for models}
  \item{estimator}{estimator upon which 'fit' is based.}
  
}
\details{Use BMA and/or model selection to form predictions using the top highest probability models.}
\seealso{\code{\link{bas}}, \code{\link{fitted.bas}}, \code{\link{confint.pred.bas}}
}
\examples{
data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")

predict(hald.gprior, newdata=Hald, estimator="BPM", se.fit=TRUE, prediction=FALSE)
# same as fitted
fitted(hald.gprior,estimator="BPM")

# default is BMA and estimation of mean vector
hald.bma = predict(hald.gprior, top=5, se.fit=TRUE)  
confint(hald.bma)

hald.bma = predict(hald.gprior, newdata=Hald[1:2,],top=5,
                    predict=TRUE, se.fit=TRUE, estimator="BMA") 
confint(hald.bma)
}
\author{Merlise Clyde}
\keyword{regression}

