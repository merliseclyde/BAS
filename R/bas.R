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




#' Bayesian Adaptive Sampling Without Replacement for Variable Selection in
#' Linear Models
#'
#' Sample without replacement from a posterior distribution on models
#'
#' BAS provides several search algorithms to find high probability models for
#' use in Bayesian Model Averaging or Bayesian model selection. For p less than
#' 20-25, BAS can enumerate all models depending on memory availability, for
#' larger p, BAS samples without replacement using random or deterministic
#' sampling. The Bayesian Adaptive Sampling algorithm of Clyde, Ghosh, Littman
#' (2010) samples models without replacement using the initial sampling
#' probabilities, and will optionally update the sampling probabilities every
#' "update" models using the estimated marginal inclusion probabilties. BAS
#' uses different methods to obtain the \code{initprobs}, which may impact the
#' results in high-dimensional problems. The deterinistic sampler provides a
#' list of the top models in order of an approximation of independence using
#' the provided \code{initprobs}.  This may be effective after running the
#' other algorithms to identify high probability models and works well if the
#' correlations of variables are small to modest.  The priors on coefficients
#' include Zellner's g-prior, the Hyper-g prior (Liang et al 2008, the
#' Zellner-Siow Cauchy prior, Empirical Bayes (local and gobal) g-priors.  AIC
#' and BIC are also included.
#'
#' @aliases bas bas.lm
#' @param formula linear model formula for the full model with all predictors,
#' Y ~ X.  All code assumes that an intercept will be included in each model
#' and that the X's will be centered.
#' @param data data frame
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#' process. Should be NULL or a numeric vector. If non-NULL, Bayes estimates
#' are obtained assuming that Y ~ N(Xb, sigma^2 diag(1/weights)).
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. The default is "na.omit".
#' @param n.models number of models to sample either without replacement
#' (method="BAS" or "MCMC+BAS") or with replacement (method="MCMC"). If NULL,
#' BAS with method="BAS" will try to enumerate all 2^p models. If enumeration
#' is not possible (memory or time) then a value should be supplied which
#' controls the number of sampled models using 'n.models'.  With method="MCMC",
#' sampling will stop once the min(n.models, MCMC.iterations) occurs so
#' MCMC.iterations be larger than n.models in order to explore the model space.
#' On exit for method= "MCMC" this is the number of unique models that have
#' been sampled with counts stored in the output as "freq""
#' @param prior prior distribution for regression coefficients.  Choices
#' include "AIC", "BIC", "g-prior", "ZS-null", "ZS-full", "hyper-g",
#' "hyper-g-laplace", "hyper-g-n", "EB-local", and "EB-global"
#' @param alpha optional hyperparameter in g-prior or hyper g-prior.  For
#' Zellner's g-prior, alpha = g, for the Liang et al hyper-g or hyper-g-n
#' method, recommended choice is alpha are between (2 < alpha < 4), with alpha
#' = 3 recommended.
#' @param modelprior Family of prior distribution on the models.  Choices
#' include \code{\link{uniform}} \code{\link{Bernoulli}} or
#' \code{\link{beta.binomial}} with the default being a
#' \code{beta.binomial(1,1)}.
#' @param initprobs Vector of length p or a character string specifiny which
#' method is used to create the vector. This is used to order variables for
#' sampling all methods for potentially more efficient storage while sampling
#' and provides the initial inclusion probabilities used for sampling without
#' replacement with method="BAS".  Options for the charactier string giving the
#' method are: "Uniform" or "uniform" where each predictor variable is equally
#' likely to be sampled (equivalent to random sampling without replacement);
#' "eplogp" uses the \code{\link{eplogprob}} function to aproximate the Bayes
#' factor from p-values from the full model to find initial marginal inclusion
#' probabilitites; "marg-eplogp" uses\code{\link{eplogprob.marg}} function to
#' aproximate the Bayes factor from p-values from the full model each simple
#' linear regression.  To run a Markov Chain to provide initial estimates of
#' marginal inclusion probabilities for "BAS", use method="MCMC+BAS" below.
#' While the initprobs are not used in sampling for method="MCMC", this
#' determines the order of the variables in the lookup table and affects memory
#' allocation in large problems where enumeration is not feasible.  For
#' variables that should always be included set the corresponding initprobs to
#' 1, i.e. the intercept should be included with probability one.
#' @param method A character variable indicating which sampling method to use:
#' method="BAS" uses Bayesian Adaptive Sampling (without replacement) using the
#' sampling probabilities given in initprobs; method="MCMC" samples with
#' replacement via a MCMC algorithm that combines the birth/death random walk
#' in Hoeting et al (1997) of MC3 with a random swap move to interchange a
#' variable in the model with one currently excluded as described in Clyde,
#' Ghosh and Littman (2010); method="MCMC+BAS" runs an initial MCMC to
#' calculate marginal inclusion probabilities and then samples without
#' replacement as in BAS.  For BAS, the sampling probabilities can be updated
#' as more models are sampled. (see update below).  We recommend "MCMC+BAS" or
#' "MCMC" for high dimensional problems where enumeration is not feasible.
#' @param update number of iterations between potential updates of the sampling
#' probabilities for method "BAS". If NULL do not update, otherwise the
#' algorithm will update using the marginal inclusion probabilities as they
#' change while sampling takes place.  For large model spaces, updating is
#' recommended. If the model space will be enumerated, leave at the default.
#' @param bestmodel optional binary vector representing a model to initialize
#' the sampling. If NULL sampling starts with the null model
#' @param prob.local A future option to allow sampling of models "near" the
#' median probability model.  Not used at this time.
#' @param prob.rw For any of the MCMC methods, probability of using the
#' random-walk proposal; otherwise use a random "flip" move to propose a new
#' model.
#' @param MCMC.iterations Number of iterations for the MCMC sampler; the
#' default is n.models*10 if not set by the user.
#' @param lambda Parameter in the AMCMC algorithm (depracated).
#' @param delta truncation parameter to prevent sampling probabilities to
#' degenerate to 0 or 1 prior to enumeration for sampling without replacement.
#' @param thin For "MCMC", thin the MCMC every "thin" iterations; default is no
#' thinning.
#' @param renormalize For MCMC sampling, should posterior probabilities be
#' based on renormalizing the marginal likelihoods times prior probabilities
#' (TRUE) or frequencies from MCMC.  The latter are unbiased in long runs,
#' while the former may have less variability.  May be compared via the
#' diagnostic plot function.
#' @return \code{bas} returns an object of class \code{bas}
#'
#' An object of class \code{BAS} is a list containing at least the following
#' components:
#'
#' \item{postprob}{the posterior probabilities of the models selected}
#' \item{priorprobs}{the prior probabilities of the models selected}
#' \item{namesx}{the names of the variables} \item{R2}{R2 values for the
#' models} \item{logmarg}{values of the log of the marginal likelihood for the
#' models} \item{n.vars}{total number of independent variables in the full
#' model, including the intercept} \item{size}{the number of independent
#' variables in each of the models, includes the intercept} \item{which}{a list
#' of lists with one list per model with variables that are included in the
#' model} \item{probne0}{the posterior probability that each variable is
#' non-zero computed using the renormalized marginal likelihoods of sampled
#' models.  This may be biased if the number of sampled models is much smaller
#' than the total number of models. Unbiased estimates may be obtained using
#' method "MCMC".} \item{mle}{list of lists with one list per model giving the
#' MLE (OLS) estimate of each (nonzero) coefficient for each model. NOTE: The
#' intercept is the mean of Y as each column of X has been centered by
#' subtracting its mean.}
#'  \item{mle.se}{list of lists with one list per model
#' giving the MLE (OLS) standard error of each coefficient for each model}
#' \item{prior}{the name of the prior that created the BMA object}
#' \item{alpha}{value of hyperparameter in prior used to create the BMA
#' object.}
#'  \item{modelprior}{the prior distribution on models that created the
#' BMA object}
#' \item{Y}{response}
#' \item{X}{matrix of predictors}
#' \item{mean.x}{vector of means for each column of X (used in
#' \code{\link{predict.bas}})}
#'
#' The function \code{\link{summary.bas}}, is used to print a summary of the
#' results. The function \code{\link{plot.bas}} is used to plot posterior
#' distributions for the coefficients and \code{\link{image.bas}} provides an
#' image of the distribution over models.  Posterior summaries of coefficients
#' can be extracted using \code{\link{coefficients.bas}}.  Fitted values and
#' predictions can be obtained using the S3 functions \code{\link{fitted.bas}}
#' and \code{\link{predict.bas}}.  BAS objects may be updated to use a
#' different prior (without rerunning the sampler) using the function
#' \code{\link{update.bas}}.
#' @author Merlise Clyde (\email{clyde@@stat.duke.edu}) and Michael Littman
#' @seealso \code{\link{summary.bas}}, \code{\link{coefficients.bas}},
#' \code{\link{print.bas}}, \code{\link{predict.bas}}, \code{\link{fitted.bas}}
#' \code{\link{plot.bas}}, \code{\link{image.bas}}, \code{\link{eplogprob}},
#' \code{\link{update.bas}}
#' @references Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive
#' Sampling for Variable Selection and Model Averaging. Journal of
#' Computational Graphics and Statistics.  20:80-101 \cr
#' \url{http://dx.doi.org/10.1198/jcgs.2010.09049}
#'
#' Clyde, M. and George, E. I. (2004) Model Uncertainty. Statist. Sci., 19,
#' 81-94. \cr \url{http://dx.doi.org/10.1214/088342304000000035}
#'
#' Clyde, M. (1999) Bayesian Model Averaging and Model Search Strategies (with
#' discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O.
#' Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.
#'
#' Hoeting, J. A., Madigan, D., Raftery, A. E. and Volinsky, C. T. (1999)
#' Bayesian model averaging: a tutorial (with discussion). Statist. Sci., 14,
#' 382-401. \cr
#' \url{http://www.stat.washington.edu/www/research/online/hoeting1999.pdf}
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2008) Mixtures
#' of g-priors for Bayesian Variable Selection. Journal of the American
#' Statistical Association.  103:410-423.  \cr
#' \url{http://dx.doi.org/10.1198/016214507000001337}
#'
#' Zellner, A. (1986) On assessing prior distributions and Bayesian regression
#' analysis with g-prior distributions. In Bayesian Inference and Decision
#' Techniques: Essays in Honor of Bruno de Finetti, pp. 233-243.
#' North-Holland/Elsevier.
#'
#' Zellner, A. and Siow, A. (1980) Posterior odds ratios for selected
#' regression hypotheses. In Bayesian Statistics: Proceedings of the First
#' International Meeting held in Valencia (Spain), pp. 585-603.
#' @keywords regression
#' @family bas methods
#' @examples
#'
#' library(MASS)
#' data(UScrime)
#' crime.bic =  bas.lm(log(y) ~ log(M) + So + log(Ed) +
#'                     log(Po1) + log(Po2) +
#'                     log(LF) + log(M.F) + log(Pop) + log(NW) +
#'                     log(U1) + log(U2) + log(GDP) + log(Ineq) +
#'                     log(Prob) + log(Time),
#'                     data=UScrime, n.models=2^15, prior="BIC",
#'                     modelprior=beta.binomial(1,1),
#'                     initprobs= "eplogp")
#'
#'
#' # use MCMC rather than enumeration
#' crime.mcmc =  bas.lm(log(y) ~ log(M) + So + log(Ed) +
#'                     log(Po1) + log(Po2) +
#'                     log(LF) + log(M.F) + log(Pop) + log(NW) +
#'                     log(U1) + log(U2) + log(GDP) + log(Ineq) +
#'                     log(Prob) + log(Time),
#'                     data=UScrime,
#'                     method="MCMC",
#'                     MCMC.iterationss=20000, prior="BIC",
#'                     modelprior=beta.binomial(1,1),
#'                     initprobs= "eplogp")
#'
#' summary(crime.bic)
#' plot(crime.bic)
#' image(crime.bic, subset=-1)
#' # more complete demo's
#' demo(BAS.hald)
#' \dontrun{demo(BAS.USCrime) }
#'
#' @rdname bas.lm
#' @family BAS methods
#' @concept BMA
#' @concept variable selection
#' @export
bas.lm = function(formula, data,  subset, weights, na.action="na.omit",
    n.models=NULL,  prior="ZS-null", alpha=NULL,
    modelprior=beta.binomial(1,1),
    initprobs="Uniform", method="BAS", update=NULL,
    bestmodel=NULL, prob.local=0.0,
    prob.rw=0.5,
    MCMC.iterations=NULL,
    lambda=NULL, delta=0.025, thin=1, renormalize=FALSE)  {


  num.updates=10
  call = match.call()

  # from lm
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())


  #data = model.frame(formula, data, na.action=na.action, weights=weights)
  n.NA = length(attr(mf, 'na.action'))

  if (n.NA > 0) {
    warning(paste("dropping ", as.character(n.NA),
                  "rows due to missing data"))
  }

  Y = model.response(mf, "numeric")
  mt <- attr(mf, "terms")
  X = model.matrix(mt, mf, contrasts)
  #X = model.matrix(formula, mf)


  Xorg = X
  namesx = dimnames(X)[[2]]
  namesx[1] = "Intercept"
  n <- dim(X)[1]

 weights=as.vector(model.weights(mf))
 if (is.null(weights)) weights = rep(1, n)

  if (length(weights) != n) stop(simpleError(paste("weights are of length ", length(weights), "not of length ", n)))

  mean.x = apply(X[,-1, drop=F], 2, weighted.mean, w=weights)
  ones = X[,1]
  X = cbind(ones, sweep(X[, -1, drop=FALSE], 2, mean.x))
  p <-  dim(X)[2]  # with intercept




  if (n <= p) {
      if (modelprior$family == "Uniform" || modelprior$family == "Bernoulli")
          warning("Uniform prior (Bernoulli)  distribution on the Model Space are not recommended for p > n; please consider using tr.beta.binomial or power.prior instead")
  }
  if (!is.numeric(initprobs)) {
      if (n <= p && initprobs == "eplogp") {
          simpleError("error: Full model is not full rank so cannot use the eplogp bound to create starting sampling probabilities, perhpas use 'marg-eplogp' for fiting marginal models\n")
      }
    initprobs = switch(initprobs,
        "eplogp" = eplogprob(lm(Y ~ X)),
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
  if (is.null(MCMC.iterations)) MCMC.iterations = as.integer(n.models*10)
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
}
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
    "BAS" = .Call(C_sampleworep_new,
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int),
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      plocal=as.numeric(prob.local),
      PACKAGE="BAS"),
    "MCMC+BAS"= .Call(C_mcmcbas,
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int),
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations),
      as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta)),
    "MCMC"= .Call(C_mcmc_new,
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int),
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations),
        as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
        as.integer(thin)),
    "AMCMC" = .Call(C_amcmc,
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int),
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      plocal=as.numeric(1.0-prob.rw), as.integer(Burnin.iterations),
      as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta)),
     "deterministic" = .Call(C_deterministic,
      Yvec, X, sqrt(weights),
      prob, modeldim,
      incint=as.integer(int),
      alpha= as.numeric(alpha),
      method=as.integer(method.num),modelprior=modelprior)
  )

  result$namesx=namesx
  result$n=length(Yvec)
  result$prior=prior
  result$modelprior=modelprior
  result$alpha=alpha
  result$probne0.RN = result$probne0
  result$postprobs.RN = result$postprobs

  if (method == "MCMC" || method == "MCMC_new" ) {
	  result$n.models = result$n.Unique
	  result$postprobs.MCMC = result$freq/sum(result$freq)

	  if (!renormalize)  {
	    result$probne0 = result$probne0.MCMC
  	  result$postprobs = result$postprobs.MCMC
	  }
  } else {
  	result$n.models=n.models
  }
  df = rep(n - 1, result$n.models)
  if (prior == "AIC" | prior == "BIC" | prior=="IC") df = df - result$size + 1
  result$df = df
  result$n.vars=p
  result$Y=Yvec
  result$X=Xorg
  result$mean.x = mean.x
  result$call=call

  result$contrasts = attr(X, "contrasts")
  result$xlevels = .getXlevels(mt, mf)
  result$terms = mt
  result$model = mf

  class(result) = c("bas")
  if (prior == "EB-global") result = EB.global(result)
  return(result)
  }
