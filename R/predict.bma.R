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


fitted.bas = function(object,  type="HPM", top=NULL, ...) {
  nmodels = length(object$which)
  X = object$X
  if (type=="HPM") {
    X = cbind(1,sweep(X[,-1], 2, object$mean.x))
    best =  min((1:nmodels)[object$logmarg == max(object$logmarg)])
    yhat  <- as.vector(X[,object$which[[best]]+1, drop=FALSE] %*% object$mle[[best]]) * object$shrinkage[[best]]
    yhat = yhat + (1 - object$shrinkage[[best]])*(object$mle[[best]])[1]
  }
  if (type == "BMA") {
   yhat = predict(object, X, top)$Ybma
}
  if (type == "MPM") {
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
}
  if (type=="BPM") {
      X = cbind(1,sweep(X[,-1], 2, object$mean.x))
      ypred = predict(object, X, top)
      sd =apply(sweep(ypred$Ypred, 2, ypred$Ybma),1, sd)
      yhat = ypred$Ypred[which.min(sd), ,drop=T]
  }
return(yhat)
}
