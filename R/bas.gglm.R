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
    		print(paste("There are", as.character(deg),
                "degenerate sampling probabilities (0 or 1); decreasing the number of models to",                 as.character(n.models)))
  	}

  	if (n.models > 2^35) stop("Dimension of model space is too big to enumerate\n  Rerun with a smaller value for n.models")
  	if (n.models > 2^20)
            print("Number of models is BIG -this may take a while")
    return(n.models)
}

.normalize.initprobs <- function (initprobs, glm.obj) {
	p <- dim(glm.obj$x)[2]
	if (!is.numeric(initprobs)) {
    		initprobs = switch(initprobs,
     			"eplogp" = eplogprob(glm.obj),
                    "uniform" = c(1.0, rep(.5, p-1)),
                    "Uniform" = c(1.0, rep(.5, p-1)),
                    )
            }
   	if (length(initprobs) == (p-1)) {
     		initprobs = c(1.0, initprobs)}
   	if (length(initprobs) != p) {
    		stop(simpleError(paste("length of initprobs is not", p)))
            }
	if (initprobs[1] < 1.0 | initprobs[1] > 1.0) initprobs[1] = 1.0
	# intercept is always included otherwise we get a segmentation
	# fault (relax later)
  	prob = as.numeric(initprobs)

	pval = summary(glm.obj)$coefficients[,4]
  	if (any(is.na(pval))) {
            print(paste("warning full model is rank deficient excluding varialble with p-values that are 0."))
            prob[is.na(pval)] = 0.0
  	}

	return(prob);
}

bas.glm = function(formula, data,  
    family = binomial(link = 'logit'),
    n.models=NULL,
    betaprior=CCH(alpha=.5, beta=nrow(data), s=0),
    modelprior=beta.binomial(1,1),
    initprobs="Uniform", 
    method="MCMC", 
    update=NULL, 
    bestmodel=NULL, 
    prob.rw=0.5,  
    MCMC.iterations=NULL,
    control = glm.control(), offset = rep(0, nobs), weights = rep(1, nobs), laplace=FALSE
                  )  {
    num.updates=10
    call = match.call()
    glm.obj = glm(formula, data, family = family, y=TRUE, x=TRUE)
	
    Y = glm.obj$y
    X = glm.obj$x
    namesx = dimnames(X)[[2]]
    namesx[1] = "Intercept" 
    p = dim(X)[2]
    nobs = dim(X)[1]
	
	
	if (is.null(offset)) {
		offset = rep(0, nobs);
	}
	
	if (is.null(weights)) {
		weights = rep(1, nobs);
	}
	
        prob <- .normalize.initprobs(initprobs, glm.obj)
	n.models <- .normalize.n.models(n.models, p, prob, method)
	modelprior <- .normalize.modelprior(modelprior,p)
  	
  	
  	#int = TRUE  # assume that an intercept is always included 
    	if (is.null(bestmodel)) {
    		bestmodel = as.integer(prob)
	}
	
  	if (is.null(update)) {
    		if (n.models == 2^(p-1))  update = n.models+1
                else (update = n.models/num.updates)
  	}


  	Yvec = as.numeric(Y)
  	modeldim = as.integer(rep(0, n.models))
  	n.models = as.integer(n.models)
    if (is.null(MCMC.iterations)) MCMC.iterations = as.integer(2*n.models)
	  
  	if (betaprior$family == "test.BF") {
  	   betaprior$hyperparameter$null.deviance =  glm(eval(call$formula), data=eval(call$data),
  	                                                 family=eval(call$family))$null.deviance
  	}
	#save(list = ls(), file = "temp.RData")
  	result = switch(method,
    		"MCMC"= .Call("glm_mcmc",
                    Y = Yvec, X = X,
                    Roffset = as.numeric(offset), Rweights = as.numeric(weights), 
                    Rprobinit = prob, Rmodeldim = modeldim,
                    modelprior= modelprior, betaprior=betaprior,
                    Rbestmodel= bestmodel,
                    plocal=as.numeric(1.0 - prob.rw),
                    BURNIN_Iterations = as.integer(MCMC.iterations),
                    family = family, Rcontrol = control, Rlaplace=as.integer(laplace),
                    PACKAGE="BAS"),
            "BAS" = .Call("glm_sampleworep",
      		Y = Yvec, X = X,
                    Roffset = as.numeric(offset), Rweights = as.numeric(weights),
                    Rprobinit = prob, Rmodeldim = modeldim,
                    modelprior = modelprior, betaprior=betaprior,
                    Rbestmodel= bestmodel,
                    plocal=as.numeric(1.0 - prob.rw),
                    family = family, Rcontrol = control,
                    Rupdate=as.integer(update), Rlaplace=as.integer(laplace),
      		PACKAGE="BAS"),
		"MCMC+BAS" = .Call("glm_mcmcbas",
      		Y = Yvec, X = X,
			Roffset = as.numeric(offset), Rweights = as.numeric(weights),
			Rprobinit = prob, Rmodeldim = modeldim,
      		modelprior = modelprior,  betaprior = betaprior,
			Rbestmodel= bestmodel,
			plocal=as.numeric(1.0 - prob.rw),
			BURNIN_Iterations = as.integer(MCMC.iterations),
			family = family, Rcontrol = control,
      		Rupdate=as.integer(update), Rlaplace=as.integer(laplace),
      		PACKAGE="BAS"),
		"deterministic" = .Call("glm_deterministic",
      		Y = Yvec, X = X,
			Roffset = as.numeric(offset), Rweights = as.numeric(weights),
			Rprobinit = prob, Rmodeldim = modeldim,
      		modelprior = modelprior,  betaprior = betaprior,
			family = family, Rcontrol = control, Rlaplace=as.integer(laplace),
      		PACKAGE="BAS")
  	)
	
  	result$namesx=namesx
  	result$n=length(Yvec)
  	result$modelprior=modelprior
    if (method == "MCMC") {
		result$n.models = result$n.Unique
            }
    else {
        result$n.models=n.models
    }

    result$R2 = .R2.glm.bas(result$deviance, result$size, call)
    result$n.vars=p
    result$Y=Yvec
    result$X=X
    result$call=call
    if (betaprior$family == "Jeffreys") result = .drop.null.bas(result)
    
    class(result) = c("basglm","bas", "bma")
    return(result) 
}

# Drop the null model from Jeffrey's prior

.drop.null.bas = function(object) {
n.models  = object$n.models


p = object$size
drop = (1:n.models)[p == 1]
logmarg = object$logmarg[-drop]
prior = object$priorprobs[-drop]

postprobs = .renormalize.postprobs(logmarg, log(prior))
which = which.matrix(object$which[-drop], object$n.var)

object$probne0 = postprobs %*% which 
object$postprobs=postprobs

method = eval(object$call$method)
if (method == "MCMC+BAS" | method == "MCMC") {
    object$freq = object$freq[-drop]
    object$probs.MCMC =  object$freq %*% which
}

object$priorprobs=prior
if (!is.null(object$sampleprobs)) object$sampleprobs = object$sampleprobs[-drop]
object$which = object$which[-drop]
object$logmarg = logmarg
object$deviance = object$deviance[-drop]
object$intercept = object$intercept[-drop]
object$size = object$size[-drop]
object$Q = object$Q[-drop]
object$R2 = object$R2[-drop]
object$mle = object$mle[-drop]
object$mle.se = object$mle.se[-drop]
object$shrinkage = object$shrinkage[-drop]
object$n.models = n.models - 1    

return(object)
}

.renormalize.postprobs = function(logmarg, logprior) {
    probs = logmarg + logprior
    probs = exp(probs - max(probs))
    probs = probs/sum(probs)
    return(probs)
}

.R2.glm.bas = function(deviance, size, call) {
    n.models = length(deviance)
    null.model = (1 : n.models)[size == 1]
    if (is.null(null.model)) {
        null.deviance =  glm(eval(call$formula), data=eval(call$data),
            family=eval(call$family))$null.deviance
    }
    else null.deviance = deviance[null.model]

    R2 = 1 - deviance/null.deviance
    return(R2)
}
