\name{confint.pred.bas}
\alias{confint.pred.bas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute Credible (Bayesian Confidence) Intervals for a BAS predict object}
\description{
Compute crebible intervals for in-sample or out of sample prediction or for the regression function
}

\usage{\method{confint}{pred.bas}(object, parm, level=.95, nsim=10000, ...)
}

\arguments{
  \item{object}{an object created by \code{\link{predict.bas}}}
  \item{parm}{character variable, "mean" or "pred".  If missing parm='pred'.}
  \item{level}{the nominal level of the (point-wise) credible interval}
  \item{nsim}{number of Monte Carlo simulations for sampling methods with BMA}
  \item{...}{optional arguments to pass on to next function call; none at this time.}
}  
\details{This constructs approximate 95 percent Highest Posterior Density intervals for 'pred.bas' objects.   If the estimator is based on model selection, the intervals use a Student t distribution using the estimate of g.  If the estimator is based on BMA, then nsim draws from the mixture of Student t distributions are obtained with the HPD interval obtained from the Monte Carlo draws.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a matrix with lower and upper level * 100 percent credible intevals for either the mean of the regression function or predicted values.   
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{Merlise A Clyde}

\seealso{\code{\link{predict.bas}}
}

\examples{
data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
hald.pred = predict(hald.gprior, estimator="BPM", predict=FALSE, se.fit=TRUE) 
confint(hald.pred, parm="mean") 
confint(hald.pred)  #default
hald.pred = predict(hald.gprior, estimator="BMA", predict=FALSE, se.fit=TRUE) 
confint(hald.pred)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")

