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
  	}

  	if (n.models > 2^30) stop("Dimension of model space is too big to enumerate\n  Rerun with a smaller value for n.models or use MCMC")
  	if (n.models > 2^25)
            warning("Number of models is BIG - this may take a while and you may run out of physical memory; you may want to consider using MCMC if your machine has limited memory.")
    return(n.models)
}




#' Bayesian Adaptive Sampling for Bayesian Model Averaging and Variable Selection in
#' Linear Models
#'
#' Sample without replacement from a posterior distribution on models
#'
#' BAS provides several algorithms to sample from posterior distributions
#' of models for
#' use in Bayesian Model Averaging or Bayesian variable selection. For p less than
#' 20-25, BAS can enumerate all models depending on memory availability.  As BAS saves all
#' models, MLEs, standard errors, log marginal likelihoods, prior and posterior and  probabilities
#' memory requirements grow linearly with M*p where M is the number of models
#' and p is the number of predictors.  For example, enumeration with p=21 with 2,097,152 takes just under
#' 2 Gigabytes on a 64 bit machine to store all summaries that would be needed for model averaging.
#' (A future version will likely include an option to not store all summaries if
#' users do not plan on using  model averaging or model selection on Best Predictive models.)
#' For larger p, BAS samples without replacement using random or deterministic
#' sampling. The Bayesian Adaptive Sampling algorithm of Clyde, Ghosh, Littman
#' (2010) samples models without replacement using the initial sampling
#' probabilities, and will optionally update the sampling probabilities every
#' "update" models using the estimated marginal inclusion probabilties. BAS
#' uses different methods to obtain the \code{initprobs}, which may impact the
#' results in high-dimensional problems. The deterministic sampler provides a
#' list of the top models in order of an approximation of independence using
#' the provided \code{initprobs}.  This may be effective after running the
#' other algorithms to identify high probability models and works well if the
#' correlations of variables are small to modest.
#' We recommend "MCMC" for
#' problems where enumeration is not feasible (memory or time constrained)
#' or even modest p if the number of
#' models sampled is not close to the number of possible models and/or there are significant
#' correlations among the predictors as the bias in estimates of inclusion
#' probabilities from "BAS" or "MSMS+BAS" may be large relative to the reduced
#' variability from using the normalized model probabilities as shown in Clyde and Ghosh, 2012.
#' Diagnostic plots with MCMC can be used to assess convergence.
#' For large problems we recommend thinning with MCMC to reduce memory requirements.
#' The priors on coefficients
#' include Zellner's g-prior, the Hyper-g prior (Liang et al 2008, the
#' Zellner-Siow Cauchy prior, Empirical Bayes (local and global) g-priors.  AIC
#' and BIC are also included, while a range of priors on the model space are available.
#'
#' @aliases bas bas.lm
#' @param formula linear model formula for the full model with all predictors,
#' Y ~ X.  All code assumes that an intercept will be included in each model
#' and that the X's will be centered.
#' @param data a data frame.  Factors will be converted to numerical vectors based on
#' the using `model.matrix`.
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
#' MCMC.iterations be significantly larger than n.models in order to explore the model space.
#' On exit for method= "MCMC" this is the number of unique models that have
#' been sampled with counts stored in the output as "freq".
#' @param prior prior distribution for regression coefficients.  Choices
#' include
#' \itemize{
#' \item "AIC"
#' \item "BIC"
#' \item "g-prior", Zellner's g prior where `g` is specified using the argument `alpha`
#' \item "JZS"  Jeffreys-Zellner-Siow prior which uses the Jeffreys
#' prior on sigma and the Zellner-Siow Cauchy prior on the coefficients.
#' The optional parameter `alpha` can be used to control
#' the squared scale of the prior, where the default is alpha=1. Setting
#' `alpha` is equal to rscale^2 in the BayesFactor package of Morey.
#' This uses QUADMATH for numerical integration of g.
#' \item  "ZS-null", a Laplace approximation to the 'JZS' prior
#' for integration of g.  alpha = 1 only. We recommend
#' using 'JZS' for accuracy and compatibility
#' with the BayesFactor package, although it is
#' slower.
#' \item "ZS-full" (to be deprecated)
#' \item "hyper-g", a mixture of g-priors where the prior on
#' g/(1+g) is a Beta(1, alpha/2) as in Liang et al (2008).  This
#' uses the Cephes library for evaluation of the marginal
#' likelihoods and may be numerically unstable for
#' large n or R2 close to 1.  Default choice of alpha is 3.
#' \item "hyper-g-laplace", Same as above but using a Laplace
#' approximation to integrate over the prior on g.
#' \item "hyper-g-n", a mixture of g-priors that where
#' u = g/n and u ~ Beta(1, alpha/2)  to provide consistency
#' when the null model is true.
#' \item "EB-local", use the MLE of g from the marginal likelihood
#' within each model
#' \item "EB-global" uses an EM algorithm to find a common or
#' global estimate of g, averaged over all models.  When it is not possible to
#' enumerate all models, the EM algorithm uses only the
#' models sampled under EB-local.
#' }
#' @param alpha optional hyperparameter in g-prior or hyper g-prior.  For
#' Zellner's g-prior, alpha = g, for the Liang et al hyper-g or hyper-g-n
#' method, recommended choice is alpha are between (2 < alpha < 4), with alpha
#' = 3 the default.  For the Zellner-Siow prior alpha = 1 by default, but can be used
#' to modify the rate parameter in the gamma prior on g,  1/g ~ G(1/2, n*alpha/2) so that
#' beta ~ C(0, sigma^2 alpha (X'X/n)^{-1}).
#' @param modelprior Family of prior distribution on the models.  Choices
#' include \code{\link{uniform}} \code{\link{Bernoulli}} or
#' \code{\link{beta.binomial}}, \code{\link{tr.beta.binomial}},
#' (with truncation) \code{\link{tr.poisson}} (a truncated Poisson), and
#' \code{\link{tr.power.prior}} (a truncated power family),
#'  with the default being a
#' \code{beta.binomial(1,1)}.  Truncated versions are useful for p > n.
#' @param initprobs Vector of length p or a character string specifying which
#' method is used to create the vector. This is used to order variables for
#' sampling all methods for potentially more efficient storage while sampling
#' and provides the initial inclusion probabilities used for sampling without
#' replacement with method="BAS".  Options for the character string giving the
#' method are: "Uniform" or "uniform" where each predictor variable is equally
#' likely to be sampled (equivalent to random sampling without replacement);
#' "eplogp" uses the \code{\link{eplogprob}} function to approximate the Bayes
#' factor from p-values from the full model to find initial marginal inclusion
#' probabilities; "marg-eplogp" uses\code{\link{eplogprob.marg}} function to
#' aproximate the Bayes factor from p-values from the full model each simple
#' linear regression.  To run a Markov Chain to provide initial estimates of
#' marginal inclusion probabilities for "BAS", use method="MCMC+BAS" below.
#' While the initprobs are not used in sampling for method="MCMC", this
#' determines the order of the variables in the lookup table and affects memory
#' allocation in large problems where enumeration is not feasible.  For
#' variables that should always be included set the corresponding initprobs to
#' 1, to override the `modelprior` or use `include.always` to force these variables
#' to always be included in the model.
#' @param include.always A formula with terms that should always be included
#' in the model with probability one.  By default this is `~ 1` meaning that the
#' intercept is always included.  This will also overide any of the values in `initprobs`
#' above by setting them to 1.
#' @param method A character variable indicating which sampling method to use:
#'\itemize{
#'\item  "deterministic" uses the "top k" algorithm described in Ghosh and Clyde (2011)
#' to sample models in order of approximate probability under conditional independence
#' using the "initprobs".  This is the most efficient algorithm for enumeration.
#'\item  "BAS" uses Bayesian Adaptive Sampling (without replacement) using the
#' sampling probabilities given in initprobs under a model of conditional independence.
#' These can be updated based on estimates of the marginal inclusion probabilities.
#' \item  "MCMC" samples with
#' replacement via a MCMC algorithm that combines the birth/death random walk
#' in Hoeting et al (1997) of MC3 with a random swap move to interchange a
#' variable in the model with one currently excluded as described in Clyde,
#' Ghosh and Littman (2010).
#' \item  "MCMC+BAS" runs an initial MCMC to
#' calculate marginal inclusion probabilities and then samples without
#' replacement as in BAS.  For BAS, the sampling probabilities can be updated
#' as more models are sampled. (see update below).
#' }
#' @param update number of iterations between potential updates of the sampling
#' probabilities for method "BAS" or "MCMC+BAS". If NULL do not update, otherwise the
#' algorithm will update using the marginal inclusion probabilities as they
#' change while sampling takes place.  For large model spaces, updating is
#' recommended. If the model space will be enumerated, leave at the default.
#' @param bestmodel optional binary vector representing a model to initialize
#' the sampling. If NULL sampling starts with the null model
#' @param prob.local A future option to allow sampling of models "near" the
#' median probability model.  Not used at this time.
#' @param prob.rw For any of the MCMC methods, probability of using the
#' random-walk Metropolis proposal; otherwise use a random "flip" move
#' to propose swap a variable that is excluded with a variable in the model.
#' @param MCMC.iterations Number of iterations for the MCMC sampler; the
#' default is n.models*10 if not set by the user.
#' @param lambda Parameter in the AMCMC algorithm (deprecated).
#' @param delta truncation parameter to prevent sampling probabilities to
#' degenerate to 0 or 1 prior to enumeration for sampling without replacement.
#' @param thin For "MCMC", thin the MCMC chain every "thin" iterations; default is no
#' thinning.  For large p, thinning can be used to significantly reduce memory
#' requirements as models and associated summaries are saved only every thin iterations.  For thin = p, the  model and associated output are recorded every p iterations,
#' similar to the Gibbs sampler in SSVS.
#' @param renormalize For MCMC sampling, should posterior probabilities be
#' based on renormalizing the marginal likelihoods times prior probabilities
#' (TRUE) or frequencies from MCMC.  The latter are unbiased in long runs,
#' while the former may have less variability.  May be compared via the
#' diagnostic plot function \code{\link{diagnostics}}.
#' See details in Clyde and Ghosh (2012).
#' @param force.heredity  Logical variable to force all levels of a factor to be
#' included together and to include higher order interactions only if lower
#' order terms are included.  Currently only supported with `method='MCMC'`.
#' Default is TRUE.
#' @param pivot Logical variable to allow pivoting of columns when obtaining the
#' OLS estimates of a model so that models that are not full rank can be fit.
#' Currently coefficients that are not estimable are set to zero.  Use caution with
#' interpreting BMA estimates of parameters.  (Experimental).
#'
#' @return \code{bas} returns an object of class \code{bas}
#'
#' An object of class \code{BAS} is a list containing at least the following
#' components:
#'
#' \item{postprob}{the posterior probabilities of the models selected}
#' \item{priorprobs}{the prior probabilities of the models selected}
#' \item{namesx}{the names of the variables}
#' \item{R2}{R2 values for the
#' models}
#' \item{logmarg}{values of the log of the marginal likelihood for the
#' models.  This is equivalent to the log Bayes Factor for comparing
#' each model to a base model with intercept only.}
#' \item{n.vars}{total number of independent variables in the full
#' model, including the intercept}
#' \item{size}{the number of independent
#' variables in each of the models, includes the intercept}
#'  \item{rank}{the rank of the design matrix; if `pivot = FALSE`, this is the same as size
#'  as no checking of rank is conducted.}
#' \item{which}{a list
#' of lists with one list per model with variables that are included in the
#' model}
#' \item{probne0}{the posterior probability that each variable is
#' non-zero computed using the renormalized marginal likelihoods of sampled
#' models.  This may be biased if the number of sampled models is much smaller
#' than the total number of models. Unbiased estimates may be obtained using
#' method "MCMC".}
#' \item{mle}{list of lists with one list per model giving the
#' MLE (OLS) estimate of each (nonzero) coefficient for each model. NOTE: The
#' intercept is the mean of Y as each column of X has been centered by
#' subtracting its mean.}
#' \item{mle.se}{list of lists with one list per model
#' giving the MLE (OLS) standard error of each coefficient for each model}
#' \item{prior}{the name of the prior that created the BMA object}
#' \item{alpha}{value of hyperparameter in coefficient prior used to create the BMA
#' object. }
#' \item{modelprior}{the prior distribution on models that created the
#' BMA object}
#' \item{Y}{response}
#' \item{X}{matrix of predictors}
#' \item{mean.x}{vector of means for each column of X (used in
#' \code{\link{predict.bas}})}
#' \item{include.always}{indices of variables that are forced into the model}
#'
#' The function \code{\link{summary.bas}}, is used to print a summary of the
#' results. The function \code{\link{plot.bas}} is used to plot posterior
#' distributions for the coefficients and \code{\link{image.bas}} provides an
#' image of the distribution over models.  Posterior summaries of coefficients
#' can be extracted using \code{\link{coefficients.bas}}.  Fitted values and
#' predictions can be obtained using the S3 functions \code{\link{fitted.bas}}
#' and \code{\link{predict.bas}}.  BAS objects may be updated to use a
#' different prior (without rerunning the sampler) using the function
#' \code{\link{update.bas}}. For MCMC sampling \code{\link{diagnostics}} can be used
#' to assess whether the MCMC has run long enough so that the posterior probabilities
#' are stable. For more details see the associated demos and vignette.
#' @author Merlise Clyde (\email{clyde@@duke.edu}) and Michael Littman
#' @seealso \code{\link{summary.bas}}, \code{\link{coefficients.bas}},
#' \code{\link{print.bas}}, \code{\link{predict.bas}}, \code{\link{fitted.bas}}
#' \code{\link{plot.bas}}, \code{\link{image.bas}}, \code{\link{eplogprob}},
#' \code{\link{update.bas}}
#' @references Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive
#' Sampling for Variable Selection and Model Averaging. Journal of
#' Computational Graphics and Statistics.  20:80-101 \cr
#' \url{http://dx.doi.org/10.1198/jcgs.2010.09049}
#'
#' Clyde, M. and Ghosh. J. (2012) Finite population estimators in stochastic search variable selection.
#' Biometrika, 99 (4), 981-988. \url{http://dx.doi.org/10.1093/biomet/ass040}
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
#' \url{http://dx.doi.org/10.1214/ss/1009212519}
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
#'
#' Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., \& Iverson, G.
#' (2009). Bayesian t-tests for accepting and rejecting the null hypothesis.
#' Psychonomic Bulletin & Review, 16, 225-237
#'
#' Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012)
#' Default Bayes Factors for ANOVA Designs. Journal of Mathematical Psychology.
#' 56.  p. 356-374.
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
#'                     initprobs= "eplogp", pivot=FALSE)
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
#'                     MCMC.iterations=20000, prior="BIC",
#'                     modelprior=beta.binomial(1,1),
#'                     initprobs= "eplogp", pivot=FALSE)
#'
#' summary(crime.bic)
#' plot(crime.bic)
#' image(crime.bic, subset=-1)
#'
#' # example with two-way interactions and hierarchical constraints
#' data(ToothGrowth)
#' ToothGrowth$dose = factor(ToothGrowth$dose)
#' levels(ToothGrowth$dose) = c("Low", "Medium", "High")
#' TG.bas = bas.lm(len ~ supp*dose, data=ToothGrowth,
#'                  modelprior=uniform(), method='BAS',
#'                  force.heredity=TRUE)
#' summary(TG.bas)
#' image(TG.bas)

#' # more complete demo's
#' demo(BAS.hald)
#' \dontrun{demo(BAS.USCrime) }
#'
#' @rdname bas.lm
#' @keywords regression
#' @family BAS methods
#' @concept BMA
#' @concept variable selection
#' @export
bas.lm = function(formula, data,  subset, weights, na.action="na.omit",
    n.models=NULL,  prior="ZS-null", alpha=NULL,
    modelprior=beta.binomial(1,1),
    initprobs="Uniform", include.always = ~ 1,
    method="BAS", update=NULL,
    bestmodel=NULL, prob.local=0.0,
    prob.rw=0.5,
    MCMC.iterations=NULL,
    lambda=NULL, delta=0.025, thin=1, renormalize=FALSE,
    force.heredity=TRUE,
    pivot=FALSE)  {


  num.updates=10
  call = match.call()
  priormethods =  c("g-prior", "hyper-g", "hyper-g-laplace", "hyper-g-n",
                    "AIC", "BIC", "ZS-null", "ZS-full",
                    "EB-local", "EB-global", "JZS")

  if (!(prior %in% priormethods)) {
    stop(paste("prior ", prior, "is not one of ",
               paste(priormethods, collapse=", ")))
  }
  # from lm
  mfall <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mfall), 0L)
  mf <- mfall[c(1L, m)]
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
        "Uniform"= c(1.0, rep(.5, p-1))
      )
  }
   if (length(initprobs) == (p-1))
       initprobs = c(1.0, initprobs)
  keep = 1
  # set up variables to always include
  if ("include.always" %in% names(mfall)) {
    minc <- match(c("include.always", "data", "subset"),  names(mfall), 0L)
    mfinc <- mfall[c(1L, minc)]
    mfinc$drop.unused.levels <- TRUE
    names(mfinc)[2] = "formula"
    mfinc[[1L]] <- quote(stats::model.frame)
    mfinc <- eval(mfinc, parent.frame())
    mtinc <- attr(mfinc, "terms")
    X.always = model.matrix(mtinc, mfinc, contrasts)

    keep = c(1L, match(colnames(X.always)[-1], colnames(X)))
    initprobs[keep] = 1.0

    if (ncol(X.always) == ncol(X)) {
      # just one model with all variables forced in
      # use method='BAS" as deterministic and MCMC fail in this context
      method='BAS'
    }
  }

  if (is.null(n.models)) n.models = min(2^p, 2^19)
  if (is.null(MCMC.iterations)) MCMC.iterations = as.integer(n.models*10)
  Burnin.iterations = as.integer(MCMC.iterations)

  if (is.null(lambda)) lambda=1.0




  int = TRUE  # assume that an intercept is always included

if (prior == "ZS-full") .Deprecated("prior='JZS'",
  msg="The Zellner-Siow full prior (Liang et al 2008)  will be deprecated in the next version of the
  package. Recommended alternative is the Jeffreys-Zellner-Siow prior 'JZS'")

# if (prior == "ZS-null") warning("We recommend using the implementation using the Jeffreys-Zellner-Siow prior (prior='JZS') which uses numerical integration rahter than the Laplace approximation")

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
      "JZS" = 9
    )

  if (is.null(alpha)) {
    alpha = switch(prior,
        "g-prior"=n,
        "hyper-g"=3,
        "EB-local"=2,
        "BIC"=n,
        "ZS-null"=1,
        "ZS-full"=n,
        "hyper-g-laplace"=3,
        "AIC"=0,
        "EB-global"=2,
        "hyper-g-n"=3,
        "JZS"= 1,
        NULL
        )
}

  if (is.null(alpha)) alpha=0.0

  parents = matrix(1,1,1)
  if (method =="MCMC+BAS" | method =="deterministic") force.heredity = FALSE # does not work with updating the tree
  if (force.heredity) {
   parents = make.parents.of.interactions(mf, data)

  # check to see if really necessary
   if (sum(parents) == nrow(parents)) {
     parents = matrix(1,1,1)
     force.heredity = FALSE}
  }

  if (is.null(bestmodel)) {
#    bestmodel = as.integer(initprobs)
    bestmodel = c(1, rep(0, p-1))
  }
    bestmodel[keep] = 1
  if (force.heredity) {
    update=NULL  # do not update tree  FIXME LATER
    if (prob.heredity(bestmodel, parents) == 0) {
      warning("bestmodel violates heredity conditions; resetting to null model")
      bestmodel = c(1, rep(0, p-1))
      }
#    initprobs=c(1, seq(.95, .55, length=(p-1) ))
  }

  prob <- .normalize.initprobs.lm(initprobs, p)
  n.models <- .normalize.n.models(n.models, p, prob, method)
  #  print(n.models)
  modelprior <- .normalize.modelprior(modelprior,p)

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
      Rparents = parents,
      Rpivot=pivot,
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
      as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
      Rparents=parents),
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
        as.integer(thin),
        Rparents=parents,
        Rpivot=pivot),
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
      method=as.integer(method.num),modelprior=modelprior, Rpivot=pivot)
  )
  result$rank_deficient = FALSE
  if (any(is.na(result$logmarg))) {
    warning("log marginals and posterior probabilities contain NA's.  Consider re-running with the option `pivot=TRUE` if there are models that are not full rank")
    result$rank_deficient=TRUE}
  if (any(result$rank != result$size)) result$rank_deficient=TRUE
  result$n.models = length(result$postprobs)
  result$namesx=namesx
  result$n=length(Yvec)
  result$prior=prior
  result$modelprior=modelprior
  result$alpha=alpha
  result$probne0.RN = result$probne0
  result$postprobs.RN = result$postprobs
  result$include.always = keep

  if (method == "MCMC" || method == "MCMC_new" ) {
	  result$n.models = result$n.Unique
	  result$postprobs.MCMC = result$freq/sum(result$freq)

	  if (!renormalize)  {
	    result$probne0 = result$probne0.MCMC
  	  result$postprobs = result$postprobs.MCMC
	  }
  }

  df = rep(n - 1, result$n.models)
  if (prior == "AIC" | prior == "BIC" | prior=="IC") df = df - result$rank + 1
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
