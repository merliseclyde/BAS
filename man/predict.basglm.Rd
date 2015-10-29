\name{predict.basglm}
\alias{predict.basglm}
\title{Prediction Method for an object of class basglm}
\description{Predictions under model averaging from a BMA (BAS) object for GLMS
}
\usage{\method{predict}{basglm}(object, newdata, top=NULL, type=c("link", "response"), ...)}

\arguments{
  \item{object}{An object of class "basglm", created by \code{bas.glm}}
  \item{newdata}{dataframe, new matrix or vector of data for predictions. May 
    include a column for the intercept or just the predictor variables.
  If a dataframe, the variables are extracted using model.matrix using
  the call that created 'object'. } 
  \item{top}{A scalar interger M.  If supplied, subset the top M models,
  based on posterior probabilities for model predictions and BMA.}
   \item{type}{Type of predictions required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable. Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities.
   }
     \item{...}{optional extra arguments}
}
\value{a list of
  \item{Ybma}{predictions using BMA}
  \item{Ypred}{matrix of predictions under each model}
    \item{postprobs}{renormalized probabilities of the top models}
  \item{best}{index of top models included}
}
\details{Use BMA to form predictions using the top highest probability
  models. 
}

\seealso{\code{\link{bas.glm}}, \code{\link{predict.bas}}, \code{\link{fitted.bas}}
}
\examples{
\dontrun{library(MASS)
data(Pima.tr)
data(Pima.te)
Pima.bas =  bas.glm(type ~ ., data=Pima.tr, family=binomial(), a=1, s=0,
b=nrow(Pima.tr)/2, modelprior=uniform(), n.models=2^7, method="BAS")
YHPM = predict(Pima.bas, newdata=Pima.te, top=1)  # Highest Probability model
}}
\author{Merlise Clyde}
\keyword{regression}
