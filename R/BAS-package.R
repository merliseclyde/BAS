#' BAS: Bayesian Model Averaging using Bayesian Adaptive Sampling
#' 
#' Package for Bayesian Model Averaging in linear models using stochastic or
#' deterministic sampling without replacement from posterior distributions.
#' Prior distributions on coefficients are of the form of Zellner's g-prior or
#' mixtures of g-priors. Options include the Zellner-Siow Cauchy Priors, the
#' Liang et al hyper-g priors, Local and Global Empirical Bayes estimates of g,
#' and other default model selection criteria such as AIC and BIC. Sampling
#' probabilities may be updated based on the sampled models.
#' 
#' \tabular{ll}{ Package: \tab BAS\cr Depends: \tab R (>= 2.8)\cr License: \tab
#' GPL-2\cr URL: http://www.stat.duke.edu/~clyde\cr }
#' 
#' Index: \preformatted{ }
#' @docType package
#' @name BAS
#' 
#' @useDynLib BAS, .registration=TRUE, .fixes="C_"
#' @aliases BAS-package BAS
#' @author Merlise Clyde, \cr Maintainer: Merlise Clyde <clyde@@stat.duke.edu>
#' @seealso \code{\link[BAS]{bas}}
#' @references Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive
#' Sampling for Variable Selection and Model Averaging. Journal of
#' Computational Graphics and Statistics.  20:80-101 \cr
#' \url{http://dx.doi.org/10.1198/jcgs.2010.09049}
#' 
#' Clyde, M. and George, E. I. (2004) Model uncertainty. Statist. Sci., 19,
#' 81-94. \cr \url{http://dx.doi.org/10.1214/088342304000000035}
#' 
#' Clyde, M. (1999) Bayesian Model Averaging and Model Search Strategies (with
#' discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O.
#' Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.
#' 
#' Li, Y. and Clyde, M. (2015) Mixtures of g-priors in Generalized Linear
#' Models.  \url{http://arxiv.org/abs/1503.06913}
#' 
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2008) Mixtures
#' of g-priors for Bayesian Variable Selection. Journal of the American
#' Statistical Association. 103:410-423.  \cr
#' \url{http://dx.doi.org/10.1198/016214507000001337}
#' @keywords package regression
#' @family bas methods
#' @examples
#' data("Hald")
#' hald.gprior =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="g-prior")
#'
#' # more complete demos
#' demo(BAS.USCrime)
#' demo(BAS.hald)
#' 
#' @import stats
#' @import graphics
#' @import grDevices
NULL







