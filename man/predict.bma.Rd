\name{predict.bas}
\alias{predict.bas}
\alias{predict}
\title{Prediction Method for an object of class BMA}
\description{Predictions under model averaging from a BMA object
}
\usage{\method{predict}{bas}(object, newdata, top=NULL, type="link", ...)}

\arguments{
  \item{object}{An object of class BAS, created by \code{bas}}
  \item{newdata}{dataframe, new matrix or vector of data for predictions. May 
    include a column for the intercept or just the predictor variables}
  \item{top}{Use only the top M models, based on posterior probabilities}
  \item{type}{Type of predictions, identity link for Gaussian models (default) provides the linear predictor}
  \item{...}{optional extra arguments}
}
\value{a list of
  \item{Ybma}{predictions using BMA}
  \item{Ypred}{matrix of predictions under each model}
  \item{best}{index of top models included}
}
\details{Use BMA to form predictions using the top highest probability
  models.
}

\seealso{\code{\link{bas}}, \code{\link{fitted.bas}}
}
\examples{
\dontrun{data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
predict(hald.gprior, hald.gprior$X[,-1], top=5)
}}
\author{Merlise Clyde}
\keyword{regression}

