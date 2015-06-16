.normalize.initprobs.lm <- function (initprobs, lm.obj) {

    p <- dim(lm.obj$x)[2]
    if (!is.numeric(initprobs)) {
    		initprobs = switch(initprobs,
     			"eplogp" = eplogprob(lm.obj),
                    "uniform" = c(1.0, rep(.5, p-1)),
                    "Uniform" = c(1.0, rep(.5, p-1)),
                    )
            }
   	if (length(initprobs) == (p-1))
     		initprobs = c(1.0, initprobs)
   	if (length(initprobs) != p)
    		stop(simpleError(paste("length of initprobs is not", p)))

	if (initprobs[1] < 1.0 | initprobs[1] > 1.0) initprobs[1] = 1.0
	# intercept is always included otherwise we get a segmentation
	# fault (relax later)
  	prob = as.numeric(initprobs)

	pval = summary(lm.obj)$coefficients[,4]
  	if (any(is.na(pval))) {
            print(paste("warning full model is rank deficient."))
#            prob[is.na(pval)] = 0.0
  	}

	return(prob);
}

bas.lm = function(formula, data, n.models=NULL,  prior="ZS-null", alpha=NULL,
                  modelprior=uniform(),
                  initprobs="Uniform", method="BAS", update=NULL, 
                  bestmodel=NULL, bestmarg=NULL, prob.local=0.0,
                  prob.rw=0.5,  
                  Burnin.iterations=NULL,
                  MCMC.iterations=NULL, lambda=NULL, delta=0.025)  {
  num.updates=10
  call = match.call()
  lm.obj = lm(formula, data, y=TRUE, x=TRUE)
  Y = lm.obj$y
  X = lm.obj$x
  Xorg = X
  namesx = dimnames(X)[[2]]
  namesx[1] = "Intercept"
  mean.x = apply(X[,-1], 2, mean)
  ones = X[,1]
  X = cbind(ones, sweep(X[, -1], 2, mean.x))
  p = dim(X)[2]

  
  if (!is.numeric(initprobs)) {
    initprobs = switch(initprobs,
     "eplogp" = eplogprob(lm.obj),
      "uniform"= c(1.0, rep(.5, p-1)),
      "Uniform"= c(1.0, rep(.5, p-1)),
      )
  }
#   if (length(initprobs) == (p-1))
#     initprobs = c(1.0, initprobs)
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
  prob <- .normalize.initprobs.lm(initprobs, lm.obj)
  n.models <- .normalize.n.models(n.models, p, prob, method)
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
  if (is.null(alpha) &&
      (method.num == 0 || method.num == 1 || method.num  == 6)) {
    stop(simpleError(paste("Must specify a value of alpha for", prior)))
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

  #MCMC-BAS
    if (is.null(MCMC.iterations)) MCMC.iterations = as.integer(n.models/2)
    if (is.null(Burnin.iterations)) Burnin.iterations = as.integer(n.models/2)
    if (is.null(lambda)) lambda=1.0

#  sampleprobs = as.double(rep(0.0, n.models))
  result = switch(method,
    "BAS" = .Call("sampleworep",
      Yvec, X,
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
      Yvec, X,
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
    "MCMC_new"= .Call("mcmc_new",
      Yvec, X,
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
    "MCMC"= .Call("mcmc",
      Yvec, X,
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
    "AMCMC" = .Call("amcmc",
      Yvec, X,
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
    "MAXeffect" = .Call("posisearch",
      Yvec, X,
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
    "deterministic" = .Call("deterministic",
      Yvec, X,
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

  class(result) = "bma"
  if (prior == "EB-global") result = EB.global.bma(result)
  return(result) 
}
