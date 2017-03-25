predict.basglm = function(object, newdata, se.fit=FALSE, 
                          type=c("response", "link"), top=NULL,
                          estimator="BMA", prediction=FALSE, ...) {
#    browser()
    if (estimator == "HPM") top=1
    
    pred = predict.bas(object, newdata, se.fit=se.fit, top=top,
                       estimator=estimator, prediction=prediction, ...)
    
    if (length(type) > 1) type = type[1]
    if (type == "response")  {
      model.specs = attributes(pred$fit)
      if (estimator == "BMA") {
        Ypred = apply(pred$Ypred, 1, 
                      FUN = function(x) {eval(object$family)$linkinv(x)})
        if (length(pred$postprobs) > 1) fit = as.vector(Ypred %*% pred$postprobs)
        else fit= as.vector(Ypred)
      }
      else fit = eval(object$family)$linkinv(pred$fit)
      attributes(fit) = model.specs
      pred$fit = fit
      if (se.fit) {
        se.fit = pred$se.fit
        se.pred = pred$se.pred
        se.fit <- se.fit * abs(eval(object$family)$mu.eta(fit))
        se.pred <- se.pred * abs(eval(object$family)$mu.eta(fit))
        pred$se.fit = se.fit
        pred$se.pred = se.pred
      }
    }

    return(pred)       
}

    
    
predict.bas = function(object, newdata, se.fit=FALSE, type="link", 
                       top=NULL,  estimator="BMA", prediction=FALSE,   ...) {
  if (!(estimator %in% c("BMA", "HPM", "MPM", "BPM"))) {
    stop("Estimator must be one of 'BMA', 'BPM', 'HPM', or 'MPM'.")
  }

  tt = terms(object)
  
  if (missing(newdata) || is.null(newdata))  {
    newdata= object$X
    insample=TRUE
    }
  else{
    if (is.data.frame(newdata)) {
    #  newdata = model.matrix(eval(object$call$formula), newdata) 
      Terms = delete.response(tt)
      m = model.frame(Terms, newdata, na.action = na.pass,
                      xlev = object$xlevels)
      newdata <- model.matrix(Terms, m, 
                              contrasts.arg = object$contrasts)
    
      insample=FALSE
    }
    else {
      stop("use of newdata as a vector is depricated, 
       please supply newdata as a dataframe")
      # if (is.vector(newdata)) newdata=matrix(newdata, nrow=1)  
    }
  }

#  browser()
  n <- nrow(newdata)
  if (ncol(newdata) == object$n.vars) newdata=newdata[,-1, drop=FALSE]  # drop intercept
  if (ncol(newdata) != (object$n.vars -1)) stop("Dimension of newdata does not match orginal model")
  if (!is.null(object$mean.x)) newdata = sweep(newdata, 2, object$mean.x)

  df = object$df

  
  if (estimator == "MPM" ) {
      nvar = object$n.vars -1
      bestmodel<- (0:nvar)[object$probne0 > .5]
      newX = cbind(1,newdata)
      best = 1
      models <- rep(0, nvar+1)
      models[bestmodel+1] <- 1
      if (sum(models) > 1) {
          object <- bas.lm(eval(object$call$formula),
                           data=eval(object$call$data), 
                           weights=eval(object$call$weights),
                           n.models=1, alpha=object$g,
                           initprobs=object$probne0, 
                           prior=object$prior, modelprior=object$modelprior,
                           update=NULL,bestmodel=models,
                           prob.local=.0)
          best= which.max(object$postprobs)
          fit  <- as.vector(newX[,object$which[[best]]+1, drop=FALSE] %*% object$mle[[best]]) * object$shrinkage[[best]]
          fit = fit + (1 - object$shrinkage[[best]])*(object$mle[[best]])[1]
          df = df[best]
      }
      else {
        fit = rep(nrow(newX), 1) * as.numeric(object$mle[object$size == 1])}
      models=bestmodel
      attributes(fit) = list(model = models, best=best)
      
      Ybma = fit
      Ypred = NULL
      postprobs=NULL
      best=NULL
      df= object$n - 1
  }
  else {    
  if (estimator == "HPM") top=1
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
  
 
  df = df[best]
  Ybma <- t(Ypred) %*% postprobs
  fit = Ybma
  if (estimator == "HPM") {
    models = unlist(object$which[best])
    attributes(fit) = list(model = models, best=best)  
  }
  if (estimator=="BPM") {
    dis =apply(sweep(Ypred, 2, Ybma),1, sd)
    bestBPM = which.min(dis)
    fit = Ypred[bestBPM, ]
    models = unlist(object$which[best[bestBPM]])
    best = best[bestBPM]
    df = df[best]
    attributes(fit) = list(model = models,
                            best = best)
   }
  
}
  #browser()
  se=list(se.fit=NULL, se.pred=NULL, 
          se.bma.fit=NULL, se.bma.pred=NULL)

  if (se.fit)  {
       if (estimator != "BMA") {
         se = .se.fit(fit, newdata, object, prediction, insample)   }
       else   {
         se = .se.bma(Ybma, newdata, Ypred, best, object, 
                      prediction, insample) }
    
  }
   
  out = list(fit=fit, Ybma=Ybma, Ypred=Ypred, postprobs=postprobs,
             se.fit=se$se.fit, se.pred=se$se.pred, 
             se.bma.fit=se$se.bma.fit, se.bma.pred=se$se.bma.pred, 
             df=df,
             best=best, bestmodel=models, 
             prediction=prediction, estimator=estimator)
    
  class(out) = 'pred.bas'  
  return(out)
}


fitted.bas = function(object,  type="response", estimator="BMA", top=NULL, ...) {
    if (type %in% c("HPM", "MPM", "BPM", "BMA")) {
        warning(paste("type = ", type,
                      " is being deprecated, use estimator = ", type))
       estimator = type
    }
        
  nmodels = length(object$which)
  X = object$X
  if (is.null(top)) top=nmodels
  if (estimator=="HPM") {
   yhat = predict(object, top=1, estimator="HPM", predict=FALSE)$fit
  }
  if (estimator == "BMA") {
   yhat = predict(object, top=top, estimator="BMA", predict=FALSE)$fit
  }
  if (estimator == "MPM") {
    yhat = predict(object, top=top, estimator="MPM", predict=FALSE)$fit
 }
  if (estimator=="BPM") {
      yhat = predict(object, top=top, estimator="BPM", predict=FALSE)$fit
  }
  
return(as.vector(yhat))
}

.se.fit= function(yhat, X, object, pred, insample) {

  n = object$n
  model = attr(yhat, "model")
  best = attr(yhat, "best")
  
  df = object$df[best]
  
  shrinkage= object$shrinkage[best]
  if (insample)  xiXTXxiT = hat(object$X[, model+1])  -1/n
  else {

    X = cbind(1, X[, model[-1], drop=FALSE] )
    oldX = (sweep(object$X[, -1], 2, object$mean.x))[, model[-1]]
#    browser()
    XRinv = X %*% solve(qr.R(qr(cbind(1,oldX))))
    xiXTXxiT = apply(XRinv^2, 1, sum) -1/n 
  }
  scale_fit = 1/n + object$shrinkage[best]*xiXTXxiT
  if (is.null(object$family)) family = gaussian()
  if (eval(family)$family == "gaussian") {
    ssy = var(object$Y)*(n-1)
    bayes_mse = ssy*(1 - shrinkage*object$R2[best])/df
    }
  else bayes_mse = 1    # ToDo add overdispersion
  se.fit = sqrt(bayes_mse*scale_fit)
  se.pred = sqrt(bayes_mse*(1 + scale_fit))
  return(list(se.fit=se.fit, se.pred=se.pred, residual.scale=sqrt(bayes_mse)))
}

.se.bma = function(fit, Xnew, Ypred, best, object, pred, insample){

n = object$n

df = object$df[best]


shrinkage= object$shrinkage[best]
if (insample) {
  xiXTXxiT =  sapply(object$which[best], 
                     FUN=function(model, X) {
                       n = nrow(X)
                       hat(X[, model[-1]+1]) -1/n}, 
                       object$X)
}
else {
  Xnew = cbind(1,Xnew)
  Xold = cbind(1,sweep(object$X[,-1], 2, object$mean.x))
  xiXTXxiT =  sapply(object$which[best], 
                     FUN=function(model, Xnew, Xold) {
                        Xnew =  Xnew[, model+1]
                        oldX = Xold[, model+1]
                        n = nrow(Xold)                
                        XRinv = Xnew %*% solve(qr.R(qr(oldX)))
                        xiXTXxiT = apply(XRinv^2, 1, sum) -1/n 
                        }, 
                     Xnew, Xold)
  
  
}

ssy = var(object$Y)*(n-1)
bayes_mse = ssy*(1 - shrinkage*object$R2[best])/df


if (is.vector(xiXTXxiT))  xiXTXxiT = matrix(xiXTXxiT, nrow=1)
  
scale_fit = 1/n + sweep(xiXTXxiT, 2, shrinkage, FUN="*")
var.fit = sweep(scale_fit, 2, bayes_mse, FUN="*")
var.pred = sweep((1 + scale_fit), 2, bayes_mse, FUN="*")


postprobs = object$postprobs[best]

# expected variance
evar.fit = as.vector(var.fit %*% postprobs)
evar.pred = as.vector(var.pred %*% postprobs)
# variance of expectations
var.efit = as.vector(postprobs %*% (sweep(Ypred, 2, fit))^2 )

se.fit = sqrt(evar.fit + var.efit)
se.pred = sqrt(evar.pred + var.efit)
  

return(list(se.bma.fit = se.fit, se.bma.pred=se.pred, 
            se.fit=t(sqrt(var.fit)), se.pred=t(sqrt(var.pred)), 
            residual.scale=sqrt(bayes_mse)))
}