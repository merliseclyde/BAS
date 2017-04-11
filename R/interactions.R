make.parents.of.interactions =
  function (mf, data) 
  { 
    termnamesX = attr(terms(mf, data=data), "term.labels") 
    p = length(termnamesX)
    interactions = grep(":", termnamesX)
    parents = diag(p)
    colnames(parents) = termnamesX
    rownames(parents) = termnamesX
    for (j in interactions) {
      term = interactions[j]
      main = unlist(strsplit(termnamesX[j], ":", 
                                      fixed = TRUE))
        parents.of.term = main
       for (i in 2:length(main)) {
        parents.of.term = c(parents.of.term, 
                            combn(main, i, FUN=paste0, collapse=":"))
              }
      parents[j, parents.of.term] = 1
    }
    
    X = model.matrix(mf, data)
    loc = attr(X, "assign")[-1] #drop intercept

    parents = parents[loc,]
    row.names(parents) = colnames(X)[-1]
    return(list(X = X, parents=parents))
  }


# model.matrix(mf, data) 
# attr( , "assign") has where terms are located

# mp = .make.parents.of.interactions(mf, df)



prob.heredity = function(model, parents, prob=.5) {
  p = length(model) 
  got.parents =  apply(parents, 1, 
           FUN=function(x){
           all(as.logical(model[as.logical(x)]))}
  )
  model.prob=0
#  browser()
  if ( all(model == got.parents)) {
    model.prob = exp(
      sum(model* log(prob) + (1 - model)*log(1.0 - prob)))
  }
  return(model.prob)
}



force_heredity.bas = function(object, prior.prob=.5) {
    mf <- object$call
    parents = make.parents.of.interactions(eval(mf$formula), eval(mf$data))$parents
    which = which.matrix(object$which, object$n.vars)
    priorprobs = apply(which[,-1], 1,
                  FUN=function(x) {prob.hereditary(model=x, parents=parents)}
                  )
    keep = (prior != 0)
    object$n.models= sum(keep)
    object$sampleprobs = object$sampleprobs[keep]   # if method=MCMC ??? reweight
    object$which = object$which[keep]
    wts = priorprobs[keep]/object$priorprobs[keep]
    method = object$call$method
    if (!is.null(method)) { 
      if (method == "MCMC" || method == "MCMC_new" ) {
         object$freq = object$freq[keep]  
         object$postprobs.MCMC = object$freq[keep]*wts
         object$postprobs.MCMC =  object$postprobs.MCMC/sum(object$postprobs.MCMC)
        object$probne0.MCMC = as.vector(postprobs.MCMC %*% which[keep,]) 
      }}
    object$priorprobs=priorprobs[keep]/sum(priorprobs[keep])
    object$logmarg = object$logmarg[keep]
    object$shrinkage=object$shrinkage[keep]
    postprobs.RN = exp(object$logmarg - min(object$logmarg))*object$priorprobs
    object$postprobs.RN = postprobs.RN/sum(postprobs.RN)
    browser()
    object$probne0.RN = as.vector(object$postprobs.RN %*% which[keep,])

    object$postprobs = object$postprobs[keep]*wts/sum(object$postprobs[keep]*wts)
    object$probne0 = as.vector(object$postprobs %*% which[keep,])

    object$mle = object$mle[keep]
    object$mle.se = object$mle.se[keep]
    object$mse = object$mse[keep]
    object$size = object$size[keep]
    object$R2 = object$R2[keep]
    object$df = object$df[keep]

  return(object)
}


data(Hald)
par.Hald = make.parents.of.interactions(Y ~ .^2, data=Hald)
bas.hald = bas.lm(Y ~ .^2, data=Hald)
hald.models = which.matrix(bas.hald$which, n.vars=bas.hald$n.vars)

prior = apply(hald.models[,-1], 1,
              FUN=function(x) {prob.hereditary(model=x, parents=par.Hald$parents)})

prob.heredity(hald.models[1,-1], par.Hald$parents)
force_heredity.bas(bas.hald)


