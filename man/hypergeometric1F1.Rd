\name{hypergeometric1F1}
\alias{hypergeometric1F1}
\title{Confluent hypergeometric2F1 function}
\description{Compute the Confluetn Hypergeometric function:
  1F1(a,b,c,t) = Gamma(b)/(Gamma(b-a)Gamma(a)) Int_0^1 t^(b-1) (1 - t)^(b-a-1) exp(c t) dt
}
\usage{hypergeometric1F1(a,b,c, laplace=TRUE, log=TRUE)}

\arguments{
  \item{a}{arbitrary}
  \item{b}{Must be greater 0}
  \item{c}{arbitrary}
  \item{laplace}{The default is to use the Laplace approximate as the
  Cephes library  routine are sometimes is unstable for large a or s
  returning NA, Inf or negative values. }
\item{log}{if TRUE, return  log(1F1)}  
}
  \references{Cephes library hyp1f1.c
 }
\examples{hypergeometric1F1(11.14756, 0.5, 0.00175097)
}
\author{Merlise Clyde (\email{clyde@stat.duke.edu})}
\keyword{math}
