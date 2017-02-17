#coefficients = function(object, ...) {UseMethod("coefficients")}
#coefficients.default = base::coefficients

coef.bas = function(object, ...) {
  conditionalmeans = list2matrix.bas(object, "mle")
  conditionalmeans[,-1] = sweep(conditionalmeans[,-1, drop=F], 1, object$shrinkage,
                                FUN="*") 
  postmean = as.vector(object$postprobs %*% conditionalmeans)

  conditionalsd  =  list2matrix.bas(object, "mle.se")
  if (!(object$prior == "AIC" || object$prior == "BIC")) {
    conditionalsd[ , -1] = sweep(conditionalsd[ ,-1, drop=F], 1, sqrt(object$shrinkage), FUN="*") }
  
  postsd = sqrt(object$postprobs %*% conditionalsd^2   +
                object$postprobs %*% ((sweep(conditionalmeans, 2, postmean, FUN="-"))^2))
  postsd = as.vector(postsd) 
  if (is.null(object$df)) {
    df = rep(object$n, length(object$postprobs))
    if (object$prior == "BIC" | object$prior == "AIC") {df = df - object$size}
    else {df = df - 1}
  }
  else df = object$df
  
  out = list(postmean=postmean, postsd=postsd, probne0 = object$probne0,
             conditionalmeans=conditionalmeans,conditionalsd=conditionalsd,
             namesx=object$namesx, postprobs=object$postprobs,
             n.vars=object$n.vars, n.models=object$n.models, df=df)
  class(out) = 'coef.bas'
  return(out)
}

print.coef.bas = function(x, n.models=5,digits = max(3, getOption("digits") - 3), ...) {
  out = cbind(x$postmean, x$postsd, x$probne0)
  dimnames(out) = list(x$namesx, c("post mean", "post SD", "post p(B != 0)"))

  cat("\n Marginal Posterior Summaries of Coefficients: \n")
  print.default(format(out, digits = digits), print.gap = 2, 
                quote = FALSE, ...)
  invisible()
}
  
# use to be pred.summary.top ???

plot.coef.bas  = function(x, e = 1e-04, subset = 1:x$n.vars, ask=TRUE, ...) {
  plotvar = function(prob0, mixprobs, df, means, sds, name,
                     e = 1e-04, nsteps = 500, ...) {
    if (prob0 == 1) {
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
    if (prob0 < 1) {
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


