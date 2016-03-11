.extractResponse <- function(frm, dat) {
# if (length(formula) == 3){
        resp <- frm[[2]];
        fdat <- eval(resp, envir=dat);
#    }
# else {stop("Formula missing Response") }
    return(fdat)
}

.normalize.initprobs.lm <- function (initprobs, p) {


    if (!is.numeric(initprobs))
    simpleError("oops no valid method given to calculate initial probabilities")
#        {
#        initprobs = switch(initprobs,
#            "eplogp" = eplogprob(lm.obj),
#            "marg-eplogp" = eplogprob.marg(Y, X),
#            "uniform" = c(1.0, rep(.5, p-1)),
#            "Uniform" = c(1.0, rep(.5, p-1)),
#            )
#    }

    if (length(initprobs) != p)
        stop(simpleError(paste("length of initprobs is", length(initprobs), "is not same as dimensions of X", p)))

    if (initprobs[1] < 1.0 | initprobs[1] > 1.0) initprobs[1] = 1.0
	# intercept is always included otherwise we get a segmentation
	# fault (relax later)
    prob = as.numeric(initprobs)
#    if (!is.null(lm.obj)) {
#        pval = summary(lm.obj)$coefficients[,4]
#        if (any(is.na(pval))) {
#            print(paste("warning full model is rank deficient."))
#        }}
#    
    return(prob);
}

.normalize.modelprior <- function(modelprior,p) {
	if (modelprior$family == "Bernoulli") {
   		if (length(modelprior$hyper.parameters) == 1) 
      		modelprior$hyper.parameters = c(1, rep(modelprior$hyper.parameters, p-1))
    		if  (length(modelprior$hyper.parameters) == (p-1)) 
     			modelprior$hyper.parameters = c(1, modelprior$hyper.parameters)
    		if  (length(modelprior$hyper.parameters) != p)
      		stop(" Number of probabilities in Bernoulli family is not equal to the number of variables or 1")
  	}
	return(modelprior)
}

.normalize.n.models <- function(n.models, p, initprobs, method) {
    if (is.null(n.models)){
        n.models = 2^(p-1)
    }
    if (n.models > 2^(p-1)) n.models = 2^(p-1)
  	deg = sum(initprobs >= 1) + sum(initprobs <= 0)
  	if (deg > 1 & n.models == 2^(p - 1)) {
    		n.models = 2^(p - deg)
    		warning(paste("There are", as.character(deg),
                "degenerate sampling probabilities (0 or 1); decreasing the number of models to",                 as.character(n.models)))
  	}

  	if (n.models > 2^30) stop("Dimension of model space is too big to enumerate\n  Rerun with a smaller value for n.models")
  	if (n.models > 2^25)
            warning("Number of models is BIG -this may take a while")
    return(n.models)
}


bas.lm = function(formula, data, weights = NULL,
    n.models=NULL,  prior="ZS-null", alpha=NULL,
    modelprior=beta.binomial(1,1),
    initprobs="Uniform", method="BAS", update=NULL, 
    bestmodel=NULL, bestmarg=NULL, prob.local=0.0,
    prob.rw=0.5,  
    MCMC.iterations=NULL,
    lambda=NULL, delta=0.025, thin=1)  {

  num.updates=10
  call = match.call()
  if ( !is.numeric(initprobs) && initprobs == "eplogp") {
      lm.obj = lm(formula, data, y=TRUE, x=TRUE)
      Y = lm.obj$y
      X = lm.obj$x
  }
  else {
      Y = .extractResponse(formula, data)    
      X = model.matrix(formula, data)
      lm.obj = NULL
  }

  Xorg = X
  namesx = dimnames(X)[[2]]
  namesx[1] = "Intercept"
  n <- dim(X)[1]

  if (is.null(weights)) {
      weights = rep(1, n);
  }

  if (length(weights) != n) stop(simpleError(paste("weights are of length ", length(weights), "not of length ", n)))
  
  mean.x = apply(X[,-1, drop=F], 2, weighted.mean, w=weights)
  ones = X[,1]
  X = cbind(ones, sweep(X[, -1], 2, mean.x))
  p <-  dim(X)[2]  # with intercept


  
  
  if (n <= p) {
      if (modelprior$family == "Uniform" || modelprior$family == "Bernoulli")
          warning("Uniform prior (Bernoulli)  distribution on the Model Space are not recommended for p > n; please consider using tr.beta.binomial instead")
  }
  if (!is.numeric(initprobs)) {
      if (n <= p && initprobs == "eplogp") {
          simpleError("error: Full model is not full rank so cannot use the eplogp bound to create starting sampling probabilities, perhpas use 'marg-eplogp' for fiting marginal models\n")
      }
    initprobs = switch(initprobs,
        "eplogp" = eplogprob(lm.obj),
        "marg-eplogp" = eplogprob.marg(Y, X),
        "uniform"= c(1.0, rep(.5, p-1)),
        "Uniform"= c(1.0, rep(.5, p-1)),
      )
  }
   if (length(initprobs) == (p-1))
       initprobs = c(1.0, initprobs)
  
#   if (length(initprobs) != p)
#    stop(simpleError(paste("length of initprobs is not", p)))

#  pval = summary(lm.obj)$coefficients[,4]
#  if (any(is.na(pval))) {
#    print(paste("warning full model is rank deficient"))
#    initprobs[is.na(pval)] = 0.0
#  }
#
#  if (initprobs[1] < 1.0 | initprobs[1] > 1.0) initprobs[1] = 1.0
# intercept is always included otherwise we get a segmentation
# fault (relax later)

  #  prob = as.numeric(initprobs)
  #MCMC-BAS
  if (is.null(n.models)) n.models = min(2^p, 2^19)
  if (is.null(MCMC.iterations)) MCMC.iterations = as.integer(n.models*2)
  Burnin.iterations = as.integer(MCMC.iterations)
  
  if (is.null(lambda)) lambda=1.0
      
  prob <- .normalize.initprobs.lm(initprobs, p)
  n.models <- .normalize.n.models(n.models, p, prob, method)
#  print(n.models)
  modelprior <- .normalize.modelprior(modelprior,p)


#  if (modelprior$family == "Bernoulli") {
#   if (length(modelprior$hyper.parameters) == 1) 
#      modelprior$hyper.parameters = c(1, rep(modelprior$hyper.parameters, p-1))
#    if  (length(modelprior$hyper.parameters) == (p-1)) 
#     modelprior$hyper.parameters = c(1, modelprior$hyper.parameters)
#    if  (length(modelprior$hyper.parameters) != p)
#      stop(" Number of probabilities in Bernoulli family is not equal to the number of variables or 1")
#}
  
  int = TRUE  # assume that an intercept is always included 
  method.num = switch(prior,
      "g-prior"=0,
      "hyper-g"=1,
      "EB-local"=2,
      "BIC"=3,
      "ZS-null"=4,
      "ZS-full"=5,
      "hyper-g-laplace"=6,
      "AIC"=7,
      "EB-global"=2,
      "hyper-g-n"=8,
    )

  if (is.null(alpha)) {
    alpha = switch(prior, 
        "g-prior"=n,
        "hyper-g"=3,
        "EB-local"=2,
        "BIC"=n,
        "ZS-null"=n,
        "ZS-full"=n,
        "hyper-g-laplace"=3,
        "AIC"=0,
        "EB-global"=2,
        "hyper-g-n"=3,
        NULL
        )
}

  if (is.null(alpha)) alpha=0.0 
  if (is.null(bestmodel)) {
    bestmodel = as.integer(initprobs)
    bestmarg = -Inf}
  if (is.null(bestmarg)) bestmarg = 0
  if (is.null(update)) {
    if (n.models == 2^(p-1))  update = n.models+1
    else (update = n.models/num.updates)
  }

  modelindex = as.list(1:n.models)
  Yvec = as.numeric(Y)
  modeldim = as.integer(rep(0, n.models))
  n.models = as.integer(n.models)

if (method == "AMCMC") {
  warning("argument method='AMCMC' is deprecated as of version 1.1.0; please use method='MCMC' instead.", 
          call. = FALSE)
}
#  sampleprobs = as.double(rep(0.0, n.models))
  result = switch(method,
    "BAS" = .Call("sampleworep",
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int), 
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      Rbestmarg=as.numeric(bestmarg),
      plocal=as.numeric(prob.local),
      PACKAGE="BAS"), 
    "MCMC+BAS"= .Call("mcmcbas",
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int), 
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      Rbestmarg=as.numeric(bestmarg),
      plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations), 
      as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
      PACKAGE="BAS"),
    "MCMC"= .Call("mcmc_new",
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int), 
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      Rbestmarg=as.numeric(bestmarg),
      plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations), 
        as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
        as.integer(thin),
	 PACKAGE="BAS"),
    "MCMC_old"= .Call("mcmc",
        Yvec, X, sqrt(weights),
        prob, modeldim,
        incint=as.integer(int), 
        alpha= as.numeric(alpha),
        method=as.integer(method.num), modelprior=modelprior,
        update=as.integer(update),
        Rbestmodel=as.integer(bestmodel),
        Rbestmarg=as.numeric(bestmarg),
        plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations), 
        as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
        as.integer(thin),
        PACKAGE="BAS"),
    "AMCMC" = .Call("amcmc",
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int), 
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      Rbestmarg=as.numeric(bestmarg),
      plocal=as.numeric(1.0-prob.rw), as.integer(Burnin.iterations), 
      as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
      PACKAGE="BAS"),
 #    "MAXeffect" = .Call("posisearch",
 #     Yvec, X,
 #     prob, modeldim,
 #     incint=as.integer(int), 
 #     alpha= as.numeric(alpha),
 #     method=as.integer(method.num), modelprior=modelprior,
 #     update=as.integer(update),
 #     Rbestmodel=as.integer(bestmodel),
 #     Rbestmarg=as.numeric(bestmarg),
 #     plocal=as.numeric(1.0-prob.rw), as.integer(Burnin.iterations), 
 #     as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
 #     PACKAGE="BAS"),
    "deterministic" = .Call("deterministic",
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int),
      alpha= as.numeric(alpha),
      method=as.integer(method.num),modelprior=modelprior,
      PACKAGE="BAS")
  )

  result$namesx=namesx
  result$n=length(Yvec)
  result$prior=prior
  result$modelprior=modelprior
  result$alpha=alpha
  if (method == "MCMC" || method == "MCMC_new" ) {
	result$n.models = result$n.Unique
  } else {
  	result$n.models=n.models
  }
  result$n.vars=p
  result$Y=Yvec
  result$X=Xorg
  result$mean.x = mean.x
  result$call=call

  class(result) = c("bas","bma")
  if (prior == "EB-global") result = EB.global.bma(result)
  return(result) 
  }
