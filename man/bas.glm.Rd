\name{bas.glm}
\alias{bas.glm}
\title{Bayesian Adaptive Sampling Without Replacement for Variable Selection in Generalized Linear Models}
\description{Sample without replacement from a posterior distribution on models}
\usage{
bas.glm(formula, data, a, b, s, 
 family = binomial(link = "logit"), 
 n.models = NULL, modelprior = beta.binomial(1, 1), 
 initprobs = "Uniform", method = "MCMC", update = NULL,
 bestmodel = NULL, bestmarg = NULL, prob.rw = 0.5, 
 Burnin.iterations = NULL, control = glm.control(), 
 offset = rep(0, nobs), weights = rep(1, nobs))
}
\arguments{
  \item{formula}{generalized linear model formula for the full model with all
    predictors, Y ~ X.  All code assumes that an intercept will be
    included in each model.}
  \item{data}{data frame}
  \item{a}{hyperparameter a for the CH-g prior.}
  \item{b}{hyperparameter b for the CH-g prior.}
  \item{s}{hyperparameter s for the CH-g prior.}
  \item{family}{a description of the error distribution and link 
    function for exponential family;
    currently only binomial() is coded.}
  \item{n.models}{number of models to sample. If NULL, BAS will
    enumerate unless p > 25.}
  \item{modelprior}{Family of prior distribution on the models.  Choices
    include \code{\link{uniform}} \code{\link{Bernoulli}} or \code{\link{beta.binomial}.}}
  \item{initprobs}{vector of length p with the initial  inclusion
    probabilities used for sampling without replacement (the intercwept
    should be included with probability one) or a character
    string giving the method used to construct the sampling probabilities
    if "Uniform" each predictor variable is equally likely to be
    sampled (equivalent to random sampling without replacement). If
    "eplogp", use the \code{\link{eplogprob}} function to aproximate the 
    Bayes factor to find initial marginal inclusion probabilitites and
    sample without replacement using these
    inclusion probabilaties. For variables that should always be
    included set the corresponding initprobs to 1. To run a
    Markov Chain to provide initial estimates of marginal
    inclusion probabilities, use method="MCMC+BAS" below.}
  \item{method}{A character variable indicating which sampling method to
    use: method="BAS" uses Bayesian Adaptive Sampling (without
    replacement) using the sampling probabilities given in initprobs;
    method="MCMC+BAS" runs an initial MCMC to calculate marginal
    inclusion probabilities and then samples without replacement as in
    BAS; method = "deterministic" runs an deterministic sampling.  For
    BAS, the sampling probabilities can be updated as more
  models are sampled. (see  update below).  We recommend "MCMC+BAS" for
  high dimensional problems.}
  \item{update}{number of iterations between potential updates of the
    sampling probabilities. If NULL do not update, otherwise the
    algorithm will update using the marginal inclusion probabilities as
    they change while sampling takes place.  For large model spaces,
    updating is recommended. If the model space will be enumerated,
    leave at the default.}
  \item{bestmodel}{optional binary vector representing a model to
    initialize the sampling. If NULL sampling starts with the null
    model}
  \item{bestmarg}{optional value for the log marginal associated with
    the bestmodel}
  \item{prob.rw}{For any of the MCMC methods, probability of using the
    random-walk proposal; otherwise use a random "flip" move to propose
    a new model.}
  \item{Burnin.iterations}{Number of iterations to discard when using any of the MCMC options}
  \item{control}{a list of parameters that control convergence in the
    fitting process.  See the documentation for
    \code{glm.control()}}
  \item{offset}{a priori known component to be included in the linear
predictor}
  \item{weights}{optional vector of weights to be used in the fitting
    process.  SHould be NULL or a numeric vector.
}
}
\details{BAS provides several search
algorithms to find high probability models for use in Bayesian Model
Averaging or Bayesian model selection. For p less than 20-25, BAS can
enumerate all models depending on memory availability, for larger p, BAS
samples without replacement using random or deterministic sampling. The
Bayesian Adaptive Sampling algorithm of Clyde, Ghosh, Littman (2010)
samples models without replacement using the initial sampling
probabilities, and will optionally update the sampling probabilities
every "update" models using the estimated marginal inclusion
probabilties. BAS uses different methods to obtain the \code{initprobs},
which may impact the results in high-dimensional problems.
The
deterinistic sampler provides a list of the top models in order of an
approximation of independence using the provided \code{initprobs}.  This
may be effective after running the other algorithms to identify high
probability models and works well if 
the correlations of variables are small to modest.  The priors on
coefficients include Zellner's g-prior, the Hyper-g prior (Liang et al
2008, the Zellner-Siow Cauchy prior, Empirical Bayes (local and gobal)
g-priors.  AIC and BIC are also included. [Need to be modified here]}
\value{
  \code{bas.glm} returns an object of class \code{BMA}

An object of class \code{BMA} is a list containing at least the following components:

\item{postprobs}{the posterior probabilities of the models selected}
\item{priorprobs}{the prior probabilities of the models selected}
  \item{logmarg}{values of the log of the marginal likelihood for the
    models}
 \item{n.vars}{total number of independent variables in the full model,
   including the intercept}
  \item{size}{the number of independent variables in each of the models,
    includes the intercept}
  \item{which}{a list of lists with one list per model with  variables
    that are included in the model} 
  \item{probne0}{the posterior probability that each variable is non-zero}
  \item{coefficients}{list of lists with one list per model giving the GLM estimate of each (nonzero) coefficient for each model.}
  \item{se}{list of lists with one list per model giving the GLM standard error of each coefficient for each model}
   \item{deviance}{the GLM deviance for each model}
  \item{modelprior}{the prior distribution on models that created the BMA object}
  \item{Q}{the Q statistic for each model}
  \item{Y}{response}
  \item{X}{matrix of predictors}
}
\references{
  Add reference here.
}
\author{Merlise Clyde (\email{clyde@stat.duke.edu}), Quanli Wang and Yingbo Li}
\note{Needed to be edited by Merlise.
}

\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{ GLM regression }
