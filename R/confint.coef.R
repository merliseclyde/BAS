confint.coef.bas = function(object, parm, level=0.95, nsim=10000, ...) {
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

    attr(ci, "Probability") <-  level
  }
  lower = paste(as.character(round(100*(1 - level)/2, 4)), " %")
  upper = paste(as.character(round(100*(1 + level)/2, 4)), " %")
  colnames(ci) = c(lower, upper)
  rownames(ci) = object$namesx[parm]
return(ci)
}

confint.pred.bas = function(object, parm, level=0.95, nsim=10000, ...) {
  
  if (missing(parm)) parm="pred"
  if (parm == "pred")  sd=object$se.pred
  else  sd = object$se.fit
  
  if (object$estimator == "BMA") {
    n.models = length(object$postprob)
    models = sample(1:n.models, size=nsim, prob= object$postprobs, replace=TRUE)
    means = object$Ypred[models,]
    df = object$df[models]
 #   browser()
    sd = sd[models,]
    npred = length(object$fit)
    pred = matrix(rt(nsim*npred, df=df), 
                   nrow=nsim, ncol=npred, byrow=FALSE)
    pred= pred*sd + means
    ci = .HPDinterval(pred, prob=level)}
  else {
 #   browser()
    df = object$df
    means = object$fit
    tq = -qt((1 - level)/2, df= df)
    ci = cbind(means - tq*sd, means + tq*sd)
    
    attr(ci, "Probability") <-  level
  }
  lower = paste(as.character(round(100*(1 - level)/2, 4)), " %")
  upper = paste(as.character(round(100*(1 + level)/2, 4)), " %")
  colnames(ci) = c(lower, upper)
  rownames(ci) = object$namesx[parm]
  return(ci)
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