predict.bas = function(object, newdata, top=NULL, ...) {
  if (is.data.frame(newdata)) {
      newdata = model.matrix(eval(object$call$formula), newdata) 
  }
  if (is.vector(newdata)) newdata=matrix(newdata, nrow=1)    
  n <- nrow(newdata)[1]
  if (ncol(newdata) == object$n.vars) newdata=newdata[,-1, drop=FALSE]  # drop intercept
  if (ncol(newdata) != (object$n.vars -1)) stop("Dimension of newdata does not match orginal model")
  if (!is.null(object$mean.x)) newdata = sweep(newdata, 2, object$mean.x)
  # need to fix for GLMS
  postprobs <- object$postprobs
  best <- order(-postprobs)
  if (!is.null(top)) best <- best[1:top]
  models <- object$which[best]
  beta <- object$ols[best]
  gg <- object$shrinkage[best]
  postprobs <- postprobs[best]
  postprobs <- postprobs/sum(postprobs)
  M <- length(postprobs)
  Ypred <- matrix(0, M, n)
  # lm case
  if (is.null(object$intercept)) {      
  for (i in 1:M) {
    beta.m <- beta[[i]]
    model.m <- models[[i]]
    Ypred[i,] <-  (newdata[,model.m[-1],drop=FALSE] %*% beta.m[-1])*gg[i]  + beta.m[1]}
}
  else {
      beta.m <- beta[[i]]
      model.m <- models[[i]]
      Ypred[i,] <-  (newdata[,model.m[-1],drop=FALSE] %*% beta.m[-1])*gg[i] + object$intercept[i]}
}
  
  Ybma <- t(Ypred) %*% postprobs
  return(list(Ybma=Ybma, Ypred=Ypred, best=best))
}


fitted.bas = function(object,  type="HPM", top=NULL, ...) {
  nmodels = length(object$which)
  X = object$X
  if (type=="HPM") {
    X = cbind(1,sweep(X[,-1], 2, object$mean.x))
    best =  min((1:nmodels)[object$logmarg == max(object$logmarg)])
    yhat  <- as.vector(X[,object$which[[best]]+1, drop=FALSE] %*% object$ols[[best]]) * object$shrinkage[[best]]
    yhat = yhat + (1 - object$shrinkage[[best]])*(object$ols[[best]])[1]
  }
  if (type == "BMA") {
   yhat = predict(object, X, top)$Ybma
}
  if (type == "MPM") {
   nvar = ncol(X) - 1
   X = cbind(1,sweep(X[,-1], 2, object$mean.x))
   bestmodel<- (0:nvar)[object$probne0 > .5]
   best = NA
   if (nvar < 32) {
      modelnum <- sapply(object$which, bin2int)
      best <- match(bin2int(bestmodel), modelnum)
    }
   if (is.na(best)) {
     model <- rep(0, nvar+1)
     model[bestmodel+1] <- 1
     object <- bas.lm(object$Y ~ object$X[,-1], n.models=1, alpha=object$g,initprobs=object$probne0, prior=object$prior, update=NULL,bestmodel=model,prob.local=.0)
     best=1
   }
   yhat  <- as.vector(X[,object$which[[best]]+1, drop=FALSE] %*% object$ols[[best]]) * object$shrinkage[[best]]
   yhat = yhat + (1 - object$shrinkage[[best]])*(object$ols[[best]])[1]
 }
return(yhat)
}
