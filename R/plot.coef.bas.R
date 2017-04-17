#' Plots the posterior distributions of coefficients derived from Bayesian
#' model averaging
#' 
#' Displays plots of the posterior distributions of the coefficients generated
#' by Bayesian model averaging over linear regression.
#' 
#' Produces plots of the posterior distributions of the coefficients under
#' model averaging.  The posterior probability that the coefficient is zero is
#' represented by a solid line at zero, with height equal to the probability.
#' The nonzero part of the distribution is scaled so that the maximum height is
#' equal to the probability that the coefficient is nonzero.
#' 
#' The parameter \code{e} specifies the range over which the distributions are
#' to be graphed by specifying the tail probabilities that dictate the range to
#' plot over.
#' 
#' @param x object of class coef.bas
#' @param e optional numeric value specifying the range over which the
#' distributions are to be graphed.
#' @param subset optional numerical vector specifying which variables to graph
#' (including the intercept)
#' @param ask Prompt for next plot
#' @param ... other parameters to be passed to \code{plot} and \code{lines}
#' @note For mixtures of g-priors, uncertainty in g is not incorporated at this
#' time, thus results are approximate
#' @author based on function \code{plot.bic} by Ian Painter in package BMA;
#' adapted for 'bas' class by Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{ \link{coef.bas}}
#' @references Hoeting, J.A., Raftery, A.E. and Madigan, D. (1996). A method
#' for simultaneous variable selection and outlier identification in linear
#' regression. Computational Statistics and Data Analysis, 22, 251-270.
#' @keywords regression
#' @examples
#' 
#' \dontrun{library(MASS)
#' data(UScrime)
#' UScrime[,-2] = log(UScrime[,-2])
#' crime.bic = bas.lm(y ~ ., data=UScrime, n.models=2^15, prior="BIC")
#' plot(coefficients(crime.bic), ask=TRUE)
#' }
#' 
#' @rdname plot.coef
#' @family bas plots
#' @method plot coef.bas
#' @export 
#' 
plot.coef.bas  = function(x, e = 1e-04, subset = 1:x$n.vars, ask=TRUE, ...) {
  plotvar = function(prob0, mixprobs, df, means, sds, name,
                     e = 1e-04, nsteps = 500, ...) {
    
    if (prob0 == 1 | length(means) == 0) {
      xlower = -0
      xupper = 0
      xmax = 1
    }
    else {
      qmin = min(qnorm(e/2, means, sds))
      qmax = max(qnorm(1 - e/2, means, sds))
      xlower = min(qmin, 0)
      xupper = max(0, qmax)
    }
    xx = seq(xlower, xupper, length.out = nsteps)
    yy = rep(0, times = length(xx))
    maxyy = 1
    if (prob0 < 1 & length(sds) > 0) {
      yy = mixprobs %*% apply(matrix(xx, ncol=1), 1,
                              FUN=function(x, d, m, s){dt(x=(x-m)/s, df=d)/s},
                              d=df, m=means, s=sds)
      maxyy = max(yy)
    }
    
    ymax = max(prob0, 1 - prob0)
    plot(c(xlower, xupper), c(0, ymax), type = "n",
         xlab = "", ylab = "", main = name, ...)
    lines(c(0, 0), c(0, prob0), lty = 1, lwd = 3, ...)
    lines(xx, (1 - prob0) * yy/maxyy, lty = 1, lwd = 1, ...)
    invisible()
  }
  
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  df = x$df
  
  for (i in subset) {
    sel = x$conditionalmeans[,i] != 0
    prob0 = 1 - x$probne0[i]   
    mixprobs = x$postprobs[sel]/(1.0 - prob0)
    means =   x$conditionalmeans[sel, i, drop=TRUE]
    sds   =   x$conditionalsd[sel, i, drop=TRUE]
    name  = x$namesx[i]
    df.sel = df[sel]
    plotvar(prob0, mixprobs, df.sel, means, sds, name, e = e, ...)
  }
  invisible()
}
