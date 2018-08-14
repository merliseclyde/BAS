
normalize.initprobs.glm <- function(initprobs, glm.obj) {
  p <- dim(glm.obj$x)[2]
  if (!is.numeric(initprobs)) {
    initprobs <- switch(initprobs,
      "eplogp" = eplogprob(glm.obj),
      "uniform" = c(1.0, rep(.5, p - 1)),
      "Uniform" = c(1.0, rep(.5, p - 1))
    )
  }
  if (length(initprobs) == (p - 1)) {
    initprobs <- c(1.0, initprobs)
  }
  if (length(initprobs) != p) {
    stop(paste("length of initprobs is not", p))
  }
  if (initprobs[1] < 1.0 | initprobs[1] > 1.0) initprobs[1] <- 1.0
  # intercept is always included otherwise we get a segmentation
  # fault (relax later)
  prob <- as.numeric(initprobs)

  pval <- summary(glm.obj)$coefficients[, 4]
  if (any(is.na(pval))) {
    warning(paste("warning full model is rank deficient; use caution when interpreting restults."))
#   prob[is.na(pval)] <- 0.0
  }

  return(prob)
}



#' Bayesian Adaptive Sampling Without Replacement for Variable Selection in
#' Generalized Linear Models
#'
#' Sample with or without replacement from a posterior distribution on GLMs
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
#' results in high-dimensional problems. The deterministic sampler provides a
#' list of the top models in order of an approximation of independence using
#' the provided \code{initprobs}.  This may be effective after running the
#' other algorithms to identify high probability models and works well if the
#' correlations of variables are small to modest.  The priors on coefficients
#' are mixtures of g-priors that provide approximations to the power prior.
#'
#' @param formula generalized linear model formula for the full model with all
#' predictors, Y ~ X.  All code assumes that an intercept will be included in
#' each model.
#' @param family a description of the error distribution and link function for
#' exponential family; currently only `binomial()` with the logistic link and
#' `poisson()` with the log link are available.
#' @param data data frame
#' @param weights optional vector of weights to be used in the fitting process.
#' May be missing in which case weights are 1.
#' @param subset subset of data used in fitting
#' @param offset a priori known component to be included in the linear
#' predictor; by default 0.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. The default is "na.omit".
#' @param n.models number of unique models to keep. If NULL, BAS will attempt
#' to enumerate unless p > 35 or method="MCMC". For any of methods using MCMC
#' algorithms that sample with replacement, sampling will stop when the number
#' of iterations exceeds the min of 'n.models' or 'MCMC.iterations' and on exit
#' 'n.models' is updated to reflect the unique number of models that have been
#' sampled.
#' @param betaprior Prior on coefficients for model coefficients (except
#' intercept).  Options include \code{\link{g.prior}}, \code{\link{CCH}},
#' \code{\link{robust}}, \code{\link{intrinsic}}, \code{\link{beta.prime}},
#' \code{\link{EB.local}}, \code{\link{AIC}}, and \code{\link{BIC}}.
#' @param modelprior Family of prior distribution on the models.  Choices
#' include \code{\link{uniform}}, \code{\link{Bernoulli}},
#' \code{\link{beta.binomial}}, truncated Beta-Binomial,
#' \code{\link{tr.beta.binomial}}, and truncated power family
#' \code{\link{tr.power.prior}}.
#' @param initprobs vector of length p with the initial inclusion probabilities
#' used for sampling without replacement (the intercept will be included with
#' probability one and does not need to be added here) or a character string
#' giving the method used to construct the sampling probabilities if "Uniform"
#' each predictor variable is equally likely to be sampled (equivalent to
#' random sampling without replacement). If "eplogp", use the
#' \code{\link{eplogprob}} function to approximate the Bayes factor using
#' p-values to find initial marginal inclusion probabilities and sample
#' without replacement using these inclusion probabilities, which may be
#' updated using estimates of the marginal inclusion probabilities. "eplogp"
#' assumes that MLEs from the full model exist; for problems where that is not
#' the case or 'p' is large, initial sampling probabilities may be obtained
#' using \code{\link{eplogprob.marg}} which fits a model to each predictor
#' separately.  To run a Markov Chain to provide initial
#' estimates of marginal inclusion probabilities, use method="MCMC+BAS" below.
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
#' method="BAS" uses Bayesian Adaptive Sampling (without replacement) using the
#' sampling probabilities given in initprobs and updates using the marginal
#' inclusion probabilities to direct the search/sample; method="MCMC" combines
#' a random walk Metropolis Hastings (as in MC3 of Raftery et al 1997) with a
#' random swap of a variable included with a variable that is currently
#' excluded (see Clyde, Ghosh, and Littman (2010) for details);
#' method="MCMC+BAS" runs an initial MCMC as above to calculate marginal
#' inclusion probabilities and then samples without replacement as in BAS;
#' method = "deterministic" runs an deterministic sampling using the initial
#' probabilities (no updating); this is recommended for fast enumeration or if a
#' model of independence is a good approximation to the joint posterior
#' distribution of the model indicators.  For BAS, the sampling probabilities
#' can be updated as more models are sampled. (see 'update' below).  We
#' recommend "MCMC+BAS" or "MCMC" for high dimensional problems.
#' @param update number of iterations between potential updates of the sampling
#' probabilities in the "BAS" method. If NULL do not update, otherwise the
#' algorithm will update using the marginal inclusion probabilities as they
#' change while sampling takes place.  For large model spaces, updating is
#' recommended. If the model space will be enumerated, leave at the default.
#' @param bestmodel optional binary vector representing a model to initialize
#' the sampling. If NULL sampling starts with the null model
#' @param prob.rw For any of the MCMC methods, probability of using the
#' random-walk proposal; otherwise use a random "flip" move to propose a new
#' model.
#' @param MCMC.iterations Number of models to sample when using any of the MCMC
#' options; should be greater than 'n.models'. By default 10*n.models.
#' @param control a list of parameters that control convergence in the fitting
#' process.  See the documentation for \code{glm.control()}
#' @param laplace logical variable for whether to use a Laplace approximate for
#' integration with respect to g to obtain the marginal likelihood.  If FALSE
#' the Cephes library is used which may be inaccurate for large n or large
#' values of the Wald Chisquared statistic.
#' @param renormalize logical variable for whether posterior probabilities
#' should be based on renormalizing marginal likelihoods times prior
#' probabilities or use Monte Carlo frequencies. Applies only to MCMC sampling.
#' @param force.heredity Logical variable to force all levels of a factor to be
#' included together and to include higher order interactions only if lower
#' order terms are included.  Currently only supported with `method='MCMC'`.
#' Default is TRUE.
#' @param bigmem Logical variable to indicate that there is access to
#' large amounts of memory (physical or virtual) for enumeration
#' with large model spaces, e.g. > 2^25.
#' @return \code{bas.glm} returns an object of class \code{basglm}
#'
#' An object of class \code{basglm} is a list containing at least the following
#' components:
#'
#' \item{postprobs}{the posterior probabilities of the models selected}
#' \item{priorprobs}{the prior probabilities of the models selected}
#' \item{logmarg}{values of the log of the marginal likelihood for the models}
#' \item{n.vars}{total number of independent variables in the full model,
#' including the intercept} \item{size}{the number of independent variables in
#' each of the models, includes the intercept} \item{which}{a list of lists
#' with one list per model with variables that are included in the model}
#' \item{probne0}{the posterior probability that each variable is non-zero}
#' \item{coefficients}{list of lists with one list per model giving the GLM
#' estimate of each (nonzero) coefficient for each model.} \item{se}{list of
#' lists with one list per model giving the GLM standard error of each
#' coefficient for each model} \item{deviance}{the GLM deviance for each model}
#' \item{modelprior}{the prior distribution on models that created the BMA
#' object} \item{Q}{the Q statistic for each model used in the marginal
#' likelihood approximation} \item{Y}{response} \item{X}{matrix of predictors}
#' \item{family}{family object from the original call} \item{betaprior}{family
#' object for prior on coefficients, including hyperparameters}
#' \item{modelprior}{family object for prior on the models}
#' \item{include.always}{indices of variables that are forced into the model}
#' @author Merlise Clyde (\email{clyde@@stat.duke.edu}), Quanli Wang and Yingbo
#' Li
#' @references Li, Y. and Clyde, M. (2015) Mixtures of g-priors in Generalized
#' Linear Models. \url{http://arxiv.org/abs/1503.06913}
#'
#' Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling for
#' Variable Selection and Model Averaging. Journal of Computational Graphics
#' and Statistics.  20:80-101 \cr
#' \url{http://dx.doi.org/10.1198/jcgs.2010.09049}
#'
#' Raftery, A.E, Madigan, D. and Hoeting, J.A. (1997) Bayesian Model Averaging
#' for Linear Regression Models. Journal of the American Statistical
#' Association.
#' @keywords GLM regression
#' @examples
#'
#' library(MASS)
#' data(Pima.tr)
#'
#'
#' # enumeration  with default method="BAS"
#' pima.cch = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
#'               method="BAS",
#'               betaprior=CCH(a=1, b=532/2, s=0), family=binomial(),
#'               modelprior=beta.binomial(1,1))
#'
#' summary(pima.cch)
#' image(pima.cch)
#'
#' pima.robust = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
#'               method="MCMC", MCMC.iterations=20000,
#'               betaprior=robust(), family=binomial(),
#'               modelprior=beta.binomial(1,1))
#'
#' pima.BIC = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
#'               method="BAS+MCMC", MCMC.iterations=1000,
#'               betaprior=bic.prior(), family=binomial(),
#'               modelprior=uniform())
#' # Poisson example
#' if(requireNamespace("glmbb", quietly=TRUE)) {
#'   data(crabs, package='glmbb')
#'   #short run for illustration
#'   crabs.bas = bas.glm(satell ~ color*spine*width + weight, data=crabs,
#'                       family=poisson(),
#'                       betaprior=EB.local(), modelprior=uniform(),
#'                       method='MCMC', n.models=2^10, MCMC.iterations=10000,
#'                       prob.rw=.95)
#' }
#' @concept BMA
#' @concept variable selection
#' @family BMA functions
#' @rdname bas.glm
#' @export
bas.glm <- function(formula, family = binomial(link = "logit"),
                    data, weights, subset, offset, na.action = "na.omit",
                    n.models = NULL,
                    betaprior = CCH(alpha = .5, beta = as.numeric(nrow(data)), s = 0),
                    modelprior = beta.binomial(1, 1),
                    initprobs = "Uniform",
                    include.always = ~1,
                    method = "MCMC",
                    update = NULL,
                    bestmodel = NULL,
                    prob.rw = 0.5,
                    MCMC.iterations = NULL,
                    control = glm.control(), laplace = FALSE, renormalize = FALSE,
                    force.heredity = TRUE,
                    bigmem = FALSE) {
  num.updates <- 10
  call <- match.call()

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if (!(family$family %in% c("binomial", "poisson"))) {
    stop(paste("family ", family$family, "not implemented"))
  }
  if (missing(data)) {
    data <- environment(formula)
  }

  # browser()
  mfall <- match.call(expand.dots = FALSE)
  m <- match(c(
    "formula", "data", "subset", "weights", "na.action",
    "etastart", "mustart", "offset"
  ), names(mfall), 0L)
  mf <- mfall[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  n.NA <- length(attr(mf, "na.action"))

  if (n.NA > 0) {
    warning(paste(
      "dropping ", as.character(n.NA),
      "rows due to missing data"
    ))
  }

  Y <- model.response(mf, type = "any")
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf, contrasts)
  # X = model.matrix(formula, mf)
  #    Y = glm.obj$y
  #    X = glm.obj$x
  namesx <- dimnames(X)[[2]]
  namesx[1] <- "Intercept"
  p <- dim(X)[2]
  nobs <- dim(X)[1]

  #   weights = as.vector(model.weights(mf))

  weights <- as.vector(model.weights(mf))
  if (is.null(weights)) {
    weights <- rep(1, nobs)
  }

  offset <- model.offset(mf)
  if (is.null(offset)) offset <- rep(0, nobs)

  #  browser()
  glm.obj <- glm(Y ~ X[, -1], family = family, weights = weights, offset = offset, y = T, x = T)

  Y <- glm.obj$y




  # set up variables to always include
  keep <- 1
  if ("include.always" %in% names(mfall)) {
    minc <- match(c("include.always", "data", "subset"), names(mfall), 0L)
    mfinc <- mfall[c(1L, minc)]
    mfinc$drop.unused.levels <- TRUE
    names(mfinc)[2] <- "formula"
    mfinc[[1L]] <- quote(stats::model.frame)
    mfinc <- eval(mfinc, parent.frame())
    mtinc <- attr(mfinc, "terms")
    X.always <- model.matrix(mtinc, mfinc, contrasts)

    keep <- c(1L, match(colnames(X.always)[-1], colnames(X)))
    prob[keep] <- 1.0
    if (ncol(X.always) == ncol(X)) {
      # just one model with all variables forced in
      # use method='BAS" as deterministic and MCMC fail in this context
      method <- "BAS"
    }
  }

  parents <- matrix(1, 1, 1)
  if (method == "deterministic" | method == "MCMC+BAS") force.heredity <- FALSE # not working yet
  if (force.heredity) {
    parents <- make.parents.of.interactions(mf, data)

    # check to see if really necessary
    if (sum(parents) == nrow(parents)) {
      parents <- matrix(1, 1, 1)
      force.heredity <- FALSE
    }
  }

  if (is.null(bestmodel)) {
    #    bestmodel = as.integer(initprobs)
    bestmodel <- c(1, rep(0, p - 1))
  }
  bestmodel[keep] <- 1
  if (force.heredity) {
    update <- NULL # do not update tree  FIXME LATER
    if (prob.heredity(bestmodel, parents) == 0) {
      warning("bestmodel violates heredity conditions; resetting to null model.  Please check  include.always and bestmodel")
      bestmodel <- c(1, rep(0, p - 1))
    }
    initprobs <- c(1, seq(.95, .55, length = (p - 1))) # keep same order
  }

  bestmodel <- as.integer(bestmodel)
  prob <- normalize.initprobs.glm(initprobs, glm.obj)
  n.models <- normalize.n.models(n.models, p, prob, method, bigmem)
  modelprior <- normalize.modelprior(modelprior, p)

  if (is.null(update)) {
    if (n.models == 2^(p - 1)) {
      update <- n.models + 1
    } else {
      (update <- n.models / num.updates)
    }
  }


  Yvec <- as.numeric(Y)
  modeldim <- as.integer(rep(0, n.models))

  if (is.null(n.models)) {
    warning("should never get here: n.models is null")
    n.models <- min(2^p, 2^19)
  }
  n.models <- as.integer(n.models)
  #print(n.models)

  if (is.null(MCMC.iterations)) {
    MCMC.iterations <- max(10000, (n.models * 10))
  }

 # MCMC.iterations = as.integer(MCMC.iterations)

  #print(MCMC.iterations)


  #  check on priors
  loglik_null <- as.numeric(-0.5 * glm(Y ~ 1,
                                       weights = weights,
                                       offset = offset,
                                       family = eval(call$family)
                                       )$null.deviance)

  betaprior$hyper.parameters$loglik_null <- loglik_null
  #  	browser()

  if (betaprior$family == "BIC" & is.null(betaprior$n)) {
    betaprior <- bic.prior(as.numeric(nobs))
  }

  if (betaprior$family == "hyper-g/n" & is.null(betaprior$n)) {
    betaprior$hyper.parameters$theta <- 1 / nobs
    betaprior$n <- nobs
  }
  if (betaprior$family == "robust" & is.null(betaprior$n)) betaprior <- robust(as.numeric(nobs))

  if (betaprior$family == "intrinsic" & is.null(betaprior$n)) {
      betaprior$hyper.parameters$n <- as.numeric(nobs)
  }

  #print(MCMC.iterations)
  result <- switch(method,
    "MCMC" = .Call(C_glm_mcmc,
      Y = Yvec, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw),
      BURNIN_Iterations = as.integer(MCMC.iterations),
      family = family, Rcontrol = control,
      Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "BAS" = .Call(C_glm_sampleworep,
      Y = Yvec, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw),
      family = family, Rcontrol = control,
      Rupdate = as.integer(update),
      Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "MCMC+BAS" = .Call(C_glm_mcmcbas,
      Y = Yvec,
      X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw),
      BURNIN_Iterations = as.integer(MCMC.iterations),
      family = family, Rcontrol = control,
      Rupdate = as.integer(update), Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "deterministic" = .Call(C_glm_deterministic,
      Y = Yvec, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      family = family,
      Rcontrol = control,
      Rlaplace = as.integer(laplace)
    )
  )

  result$namesx <- namesx
  result$n <- length(Yvec)
  result$modelprior <- modelprior
  result$probne0.RN <- result$probne0
  result$postprobs.RN <- result$postprobs
  result$family <- family
  result$betaprior <- betaprior
  result$modelprior <- modelprior

  result$n.models <- length(result$postprobs)
  result$include.always <- keep

  # 	if (method == "MCMC") result$n.models = result$n.Unique


  df <- rep(nobs - 1, result$n.models)

  if (betaprior$class == "IC") df <- df - result$size + 1
  result$df <- df
  result$R2 <- .R2.glm.bas(result$deviance, result$size, call)
  result$n.vars <- p
  result$Y <- Yvec
  result$X <- X
  result$call <- call
  result$terms <- mt
  result$contrasts <- attr(X, "contrasts")
  result$xlevels <- .getXlevels(mt, mf)
  result$model <- mf

  # drop null model
  if (betaprior$family == "Jeffreys") result <- .drop.null.bas(result)

  if (method == "MCMC") {
    result$postprobs.MCMC <- result$freq / sum(result$freq)
    if (!renormalize) {
      result$probne0 <- result$probne0.MCMC
      result$postprobs <- result$postprobs.MCMC
    }
  }

  class(result) <- c("basglm", "bas")
  return(result)
}

# Drop the null model from Jeffrey's prior

.drop.null.bas <- function(object) {
  n.models <- object$n.models


  p <- object$size
  drop <- (1:n.models)[p == 1]
  logmarg <- object$logmarg[-drop]
  prior <- object$priorprobs[-drop]

  postprobs <- .renormalize.postprobs(logmarg, log(prior))
  which <- which.matrix(object$which[-drop], object$n.var)

  object$probne0 <- postprobs %*% which
  object$postprobs <- postprobs

  method <- eval(object$call$method)
  if (method == "MCMC+BAS" | method == "MCMC") {
    object$freq <- object$freq[-drop]
    object$probne0.MCMC <- object$freq %*% which
  }

  object$priorprobs <- prior
  if (!is.null(object$sampleprobs)) object$sampleprobs <- object$sampleprobs[-drop]
  object$which <- object$which[-drop]
  object$logmarg <- logmarg
  object$deviance <- object$deviance[-drop]
  object$intercept <- object$intercept[-drop]
  object$size <- object$size[-drop]
  object$Q <- object$Q[-drop]
  object$R2 <- object$R2[-drop]
  object$mle <- object$mle[-drop]
  object$mle.se <- object$mle.se[-drop]
  object$shrinkage <- object$shrinkage[-drop]
  object$n.models <- n.models - 1
  object$df <- object$df[-drop]
  return(object)
}

.renormalize.postprobs <- function(logmarg, logprior) {
  probs <- logmarg + logprior
  probs <- exp(probs - max(probs))
  probs <- probs / sum(probs)
  return(probs)
}

.R2.glm.bas <- function(deviance, size, call) {
  n.models <- length(deviance)
  null.model <- (1:n.models)[size == 1]
  if (is.null(null.model)) {
    null.deviance <- glm(eval(call$formula),
      data = eval(call$data),
      family = eval(call$family)
    )$null.deviance
  }
  else {
    null.deviance <- deviance[null.model]
  }

  R2 <- 1 - deviance / null.deviance
  return(R2)
}
