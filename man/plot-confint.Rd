\name{plot.confint.bas}
\alias{plot.confint.bas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Bayesian Confidence Intervals}
\description{
Function takes the the output of functions that return credible intervals from BAS objects, and creates a plot of the posterior mean  with segments representing the credible interval. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{# S3 plotting method for objects of class 'confint.bas'
\method{plot}{confint.bas}(x, horizontal=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the output from \code{\link{confint.coef.bas}} or \code{\link{confint.pred.bas}} containing credible intervals and estimates.
  }
\item{horizontal}{orientation of the plot}
\item{...}{optional graphical arguments to pass on to plot}
}
\details{
This function takes the HPD intervals or credible intervals created by \code{\link{confint.coef.bas}} or \code{\link{confint.pred.bas}} from BAS objects, and creates a plot of the posterior mean  with segments representing the credible interval.  BAS tries to return HPD intervals, and under model averaging these may not be symmetric.  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot of the credible intervals.
}
\author{
Merlise A Clyde
}

\seealso{\code{\link{confint.coef.bas}}, \code{\link{confint.pred.bas}},
\code{\link{coef.bas}}, \code{\link{predict.bas}}, \code{link{bas.lm}}
}
\examples{
data(Hald)
hald.ZS = bas.lm(Y ~ ., data=Hald, prior="ZS-null", modelprior=uniform())
plot(confint(coef(hald.ZS),parm=2:5))
plot(confint(predict(hald.ZS, se.fit=TRUE), parm="mean"))
}
\keyword{regression}% use one of  RShowDoc("KEYWORDS")
\keyword{bayesian}% __ONLY ONE__ keyword per line
