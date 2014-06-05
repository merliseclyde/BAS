\name{bin2int}
\alias{bin2int}
\title{Convert binary model representation into an integer}
\description{Takes a binary string representation of a model and
  converts to an integer}

\usage{bin2int(model)}

\arguments{\item{model}{a Boolean/binary vector of length p representing
    a model}}

\value{an integer}

\details{Used in \code{fitted.bma} to determine if the median
  probability model is included in the sample.  Not meant to be used
  directly by the user.  On a 32 bit system, p must be less than or
  equal to 32.
  }


\author{Merlise Clyde \email{clyde@stat.duke.edu}}
\keyword{regression}
