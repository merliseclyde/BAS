update.bma = function(object, newprior, alpha=NULL, ...) {
  method.num = switch(newprior,
    "g-prior"=0,
    "hyper-g"=1,
    "EB-local"=2,
    "BIC"=3,
    "ZS-null"=4,
    "ZS-full"=5,
    "hyper-g-laplace"=6,
    "AIC"=7,
    "EB-global"=9,
    "hyper-g-n"=8
    )
  if (is.null(alpha) &&
      (method.num == 0 || method.num == 1 || method.num == 6)) {
    stop(paste("Must specify a value of alpha for", newprior))
  }
  
  if (is.null(alpha)) alpha=0.0
  object$alpha = alpha

  if (newprior == "EB-global") object = EB.global.bma(object)
  else {
    object$prior= newprior
    SSY = sum((object$Y - mean(object$Y))^2)
    R2Full = summary(lm(object$Y ~ object$X[,-1]))$r.squared
    logmarg=object$logmarg
    shrinkage=object$shrinkage
    tmp = .C("gexpectations_vect",  nmodels=as.integer(length(object$which)),
      p=as.integer(object$n.vars),  pmodel=as.integer(object$size),
      nobs=as.integer(object$n), R2=object$R2, alpha=as.double(alpha),
      method=as.integer(method.num), RSquareFull=as.double(R2Full), SSY=as.double(SSY),
      logmarg=logmarg, shrinkage=shrinkage,
      PACKAGE="BAS")
    object$logmarg = tmp$logmarg
    object$shrinkage=tmp$shrinkage
    object$postprobs = exp(object$logmarg - min(object$logmarg))*object$priorprobs
    object$postprobs = object$postprobs/sum(object$postprobs)
    which = which.matrix(object$which, object$n.vars)
    object$probne0 = object$postprobs %*% which
  }
  return(object)
}

