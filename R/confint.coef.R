confint.coef.bas = function(object, parm, level=0.95, nsim=10000, plot=F, ...) {
  n.models = length(object$postprob)
  if (missing(parm)) parm= 1:object$n.vars
  
  if (n.models > 1) {
    models = sample(1:n.models, size=nsim, prob= object$postprobs, replace=TRUE)
    means = object$conditionalmeans[models,parm]
    sd = object$conditionalsd[models,parm]
    df = object$df
    if (length(df) == length(object$postprobs)) df = object$df[models]
    betas = matrix(rt(nsim*length(subset), df=df), 
                   nrow=nsim, ncol=length(parm), byrow=FALSE)
    betas = betas*sd + means
    ci = .HPDinterval(betas, prob=level)}
  else {
    df = sum(object$postprobs*object$df)
    means = object$postmean[parm]
    sd = object$postsd[parm]
    tq = -qt((1 - level)/2, df= df)
    ci = cbind(means - tq*sd, means + tq*sd)
  }
  ci = cbind(ci, object$postmean[parm])
  attr(ci, "Probability") <-  level
  attr(ci, "class") = "confint.bas"
  lower = paste(as.character(round(100*(1 - level)/2, 4)), " %")
  upper = paste(as.character(round(100*(1 + level)/2, 4)), " %")
  colnames(ci) = c(lower, upper, "beta")
  rownames(ci) = object$namesx[parm]
return(ci)
}

confint.pred.bas = function(object, parm, level=0.95, nsim=10000, ...) {
  
  if (missing(parm)) parm="pred"
  if (parm == "pred")  sd=object$se.pred
  else  sd = object$se.fit
  
  if (is.null(sd)) {
    warning("object does not have fitted or prediction standard deviations")
    return()}
  if (object$estimator == "BMA") {
    n.models = length(object$postprob)
    models = sample(1:n.models, size=nsim, prob= object$postprobs, replace=TRUE)
    means = object$Ypred[models,]
    df = object$df[models]
 #  # browser()
    sd = sd[models,]
    npred = length(object$fit)
    pred = matrix(rt(nsim*npred, df=df), 
                   nrow=nsim, ncol=npred, byrow=FALSE)
    pred= pred*sd + means
    ci = .HPDinterval(pred, prob=level)
  }
  else {
 #   browser()
    df = object$df
    means = object$fit
    tq = -qt((1 - level)/2, df= df)
    ci = cbind(means - tq*sd, means + tq*sd)
  }
  
  # browser()
  ci = cbind(ci, object$fit)
  attr(ci, "Probability") <-  level
  attr(ci, "class") = "confint.bas"
  lower = paste(as.character(round(100*(1 - level)/2, 4)), " %")
  upper = paste(as.character(round(100*(1 + level)/2, 4)), " %")
  colnames(ci) = c(lower, upper, parm)
  rownames(ci) = object$namesx[parm]
  
  
  return(ci)
}

plot.confint.bas = function(ci, horizontal=FALSE, ...) {  
    namesx = rownames(ci)   
    y = ci[,3]
    x = 1:nrow(ci)
    xlim = range(x) + c(-0.5,0.2)
    ylim = range(pretty(ci))
    
    type = colnames(ci)[3]
    if (type == "beta") {
      ylab = bquote(hat(beta))
      xlab = "coefficient" }
    else {
      ylab = bquote(hat(Y))
      xlab = type
    }
    par(mar=c(4.5,5,1,1), las=1)
    if (!horizontal) {
      plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n", ...)
      axis(1, at=x, labels=namesx, tick=FALSE, ...)
      abline(h=0, lty=3, ...)
      arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05, ...)
    }
### horizontal layout:
    else { 
      plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, 
           xlab=ylab, ylab="", yaxt="n", bty="n", ...)
      axis(2, at=x, labels=namesx, tick=FALSE, ...)
      abline(v=0, lty=3, ...)
      arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05, ...)
    }
return()
}

.HPDinterval = function (obj, prob = 0.95, ...) 
{
  # from library coda but used here so that library does not have to be loaded
  obj <- as.matrix(obj)
  vals <- apply(obj, 2, sort)
  if (!is.matrix(vals)) 
    stop("obj must have nsamp > 1")
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, , drop = FALSE] - vals[init, 
                                                        , drop = FALSE], 2, which.min)
  ans <- cbind(vals[cbind(inds, 1:npar)], vals[cbind(inds + 
                                                       gap, 1:npar)])
  lower = as.character(round(100*(1 - prob)/2, 2))
  upper = as.character(round(100*(prob +1)/2, 2))
  
  dimnames(ans) <- list(colnames(obj), c(lower, upper))
  attr(ans, "Probability") <- gap/nsamp
  ans
}