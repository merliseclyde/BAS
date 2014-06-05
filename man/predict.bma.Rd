\name{predict.bma}
\alias{predict.bma}
\alias{predict}
\title{Prediction Method for an object of class BMA}
\description{Predictions under model averaging from a BMA object
}
\usage{\method{predict}{bma}(object, newdata, top=NULL, ...)}

\arguments{
  \item{object}{An object of class BMA, created by \code{bas}}
  \item{newdata}{new matrix or vector of data for predictions. May 
    include a column for the intercept or just the predictor variables}
  \item{top}{Use only the top M models, based on posterior
    probabilities}
  \item{...}{optional extra arguments}
}
\value{a list of
  \item{Ybma}{predictions using BMA}
  \item{Ypred}{matrix of predictions under each model}
  \item{best}{index of top models included}
}
\details{Use BMA to form predictions using the top highest probability
  models. Currently newdata must be in the form of a matrix or vector
  with variables in the same order as in the model matrix used to obtain
  the BMA object (see \code{object$X}).  Future versions 
  will allow newdata to be a dataframe.
}

\seealso{\code{\link{bas}}, \code{\link{fitted.bma}}
}
\examples{
\dontrun{data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
predict(hald.gprior, hald.gprior$X[,-1], top=5)
}}
\author{Merlise Clyde}
\keyword{regression}

