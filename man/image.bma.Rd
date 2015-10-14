\name{image.bma}
\alias{image.bma}
\alias{image}
\title{Images of models used in Bayesian model averaging }
\description{
  Creates an image of the models selected using \code{\link{bas}}.
}
\usage{\method{image}{bma}(x, top.models=20, intensity=TRUE, prob=TRUE, log=TRUE,
rotate=TRUE, color="rainbow", subset=NULL, offset=.75, digits=3,
vlas=2,plas=0,rlas=0, ...)
}
\arguments{
  \item{x}{An object of type 'bma' created by BAS}
  \item{top.models}{Number of the top ranked models to plot}
  \item{intensity}{Logical variable, when TRUE image intensity is
    proportional to the probability or log(probability) of the model,
    when FALSE, intensity is binary indicating just presence (light) or
    absence (dark) of a variable.}
  \item{prob}{Logical variable for whether the area in the image for
    each model should be proportional to the posterior probability (or log
    probability) of the model (TRUE) or with equal area (FALSE). }
  \item{log}{Logical variable indicating whether the intensities should
    be based on log posterior odds (TRUE) or posterior
    probabilities  (FALSE).  The log of the posterior odds is for
    comparing the each  model to the worst model in the top.models.}
  \item{rotate}{Should the image of models be rotated so that models are
    on the y-axis and variables are on the x-axis (TRUE)}
  \item{color}{The color scheme for image intensities. The value
    "rainbow" uses the rainbow palette. The value "blackandwhite"
    produces a black  and white image (greyscale image)}
  \item{subset}{indices of variables to include in plot; 1 is the
    intercept}
  \item{offset}{numeric value to add to intensity}
  \item{digits}{number of digits in posterior probabilities to keep}
  \item{vlas}{las parameter for placing variable names; see par}
  \item{plas}{las parameter for posterior probability axis}
  \item{rlas}{las parameter for model ranks}
  \item{...}{Other parameters to be passed to the \code{image} and \code{axis} functions.}}
\details{ 
  Creates an image of the model space sampled using \code{\link{bas}}.
  If a subset of the top models are plotted, then probabilities are
  renormalized over the subset.
  
  }
\author{Merlise Clyde \email{clyde@stat.duke.edu}}
\references{Clyde, M. (1999) Bayesian Model Averaging and Model Search Strategies (with discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.
}

\seealso{\code{\link{bas}}}
\examples{
require(graphics)
data("Hald")
hald.ZSprior =  bas.lm(Y~ ., data=Hald,  prior="ZS-null")
image(hald.ZSprior, subset=-1)
}
\note{Suggestion to allow area of models be proportional to posterior
  probability due to Thomas Lumley} 
\keyword{regression}

