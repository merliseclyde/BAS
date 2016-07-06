predict.basglm = function(object, newdata, top=NULL, type=c("link", "response"), ...) {
#    browser()
    pred = predict.bas(object, newdata, top)
    if (length(type) > 1) type = type[1]
    if (type == "response") {
        Ypred = apply(pred$Ypred, 1, FUN = function(x) {eval(object$call$family)$linkinv(x)})
        if (top > 1) {
            Ybma = Ypred %*% pred$postprobs}
        else Ybma = Ypred

        pred = list(Ypred=Ypred, Ybma=Ybma, postprobs=pred$postprobs, best=pred$best)
    }   
    return(pred)       
}

    
    
predict.bas = function(object, newdata, top=NULL, type="link", ...) {

  if (is.data.frame(newdata)) {
      newdata = model.matrix(eval(object$call$formula), newdata) 
  }
  if (is.vector(newdata)) newdata=matrix(newdata, nrow=1)    
  n <- nrow(newdata)
  if (ncol(newdata) == object$n.vars) newdata=newdata[,-1, drop=FALSE]  # drop intercept
  if (ncol(newdata) != (object$n.vars -1)) stop("Dimension of newdata does not match orginal model")
  if (!is.null(object$mean.x)) newdata = sweep(newdata, 2, object$mean.x)

  postprobs <- object$postprobs
  best <- order(-postprobs)
  if (!is.null(top)) best <- best[1:top]
  models <- object$which[best]
  beta <- object$mle[best]
  gg <- object$shrinkage[best]
  intercept <- object$intercept[best]  
  postprobs <- postprobs[best]
  postprobs <- postprobs/sum(postprobs)
  M <- length(postprobs)
  Ypred <- matrix(0, M, n)
                                        # lm case
  if (is.null(intercept)) {      
      for (i in 1:M) {
          beta.m <- beta[[i]]
          model.m <- models[[i]]
          Ypred[i,] <-  (newdata[,model.m[-1],drop=FALSE] %*% beta.m[-1])*gg[i]  + beta.m[1]}
  }
  else {
     for (i in 1:M) { 
      beta.m <- beta[[i]]
      model.m <- models[[i]]
      Ypred[i,] <-  (newdata[,model.m[-1],drop=FALSE] %*% beta.m[-1])*gg[i] + intercept[i]}
 }
  
  Ybma <- t(Ypred) %*% postprobs
  return(list(Ybma=Ybma, Ypred=Ypred, postprobs=postprobs, best=best))
}


fitted.bas = function(object,  type="response", estimator=NULL, top=NULL, ...) {
    if (type %in% c("HPM", "MPM", "BPM", "BMA")) {
        warning(paste("type = ", type,
                      " is being deprecated, use estimator = ", type))
        if (is.null(estimator)) estimator = type
    }
        
  nmodels = length(object$which)
  X = object$X
  if (estimator=="HPM") {
#    X = cbind(1,sweep(X[,-1], 2, object$mean.x))
#    best =  which.max(object$logmarg)
#    yhat  <- as.vector(X[,object$which[[best]]+1, drop=FALSE] %*% object$mle[[best]]) * object$shrinkage[[best]]
 #   yhat = yhat + (1 - object$shrinkage[[best]])*(object$mle[[best]])[1]
      ypred = predict(object, X, 1)
      best = ypred$best
      yhat = ypred$Ybma   # note with ome model this is the HPM
      attributes(yhat) = list(model = unlist(object$which[best]), best=best)   
  }
  if (estimator == "BMA") {
   yhat = predict(object, X, top)$Ybma
}
  if (estimator == "MPM") {
   nvar = ncol(X) - 1
   X = cbind(1,sweep(X[,-1], 2, object$mean.x))
   bestmodel<- (0:nvar)[object$probne0 > .5]
   best = NA
   model <- rep(0, nvar+1)
   model[bestmodel+1] <- 1
   if (sum(model) > 1) {
       object <- bas.lm(eval(object$call$formula), data=eval(object$call$data), n.models=1, alpha=object$g,
                        initprobs=object$probne0,
                        prior=object$prior, update=NULL,bestmodel=model,
                        prob.local=.0)
     best=1
     yhat  <- as.vector(X[,object$which[[best]]+1, drop=FALSE] %*% object$mle[[best]]) * object$shrinkage[[best]]
     yhat = yhat + (1 - object$shrinkage[[best]])*(object$mle[[best]])[1]
 }
   else { yhat = rep(nrow(X), 1) * as.numeric(object$mle[object$size == 1])}
   attributes(yhat) = list(model = bestmodel)   
}
  if (estimator=="BPM") {
      ypred = predict(object, X, top)
      dis =apply(sweep(ypred$Ypred, 2, ypred$Ybma),1, sd)
      best = which.min(dis)
      yhat = ypred$Ypred[best, ]
      attributes(yhat) = list(model = unlist(object$which[ypred$best[best]]),
                    best = ypred$best[best])
  }
return(yhat)
}
