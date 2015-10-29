\name{bas.glm}
\alias{bas.glm}
\title{Bayesian Adaptive Sampling Without Replacement for Variable Selection in Generalized Linear Models}
\description{Sample with or without replacement from a posterior distribution on GLMs}
\usage{
 bas.glm(formula, data, 
 family = binomial(link = "logit"), 
 n.models = NULL, betaprior=CCH(alpha=.5, beta=nrow(data), s=0),
 modelprior = beta.binomial(1,1), 
 initprobs = "Uniform", method = "MCMC", update = NULL,
 bestmodel = NULL, bestmarg = NULL, prob.rw = 0.5, 
 MCMC.iterations = NULL, control = glm.control(), 
 offset = rep(0, nobs), weights = rep(1, nobs), laplace=FALSE)
}
\arguments{
  \item{formula}{generalized linear model formula for the full model with all
    predictors, Y ~ X.  All code assumes that an intercept will be
    included in each model.}
  \item{data}{data frame}
  \item{family}{a description of the error distribution and link 
    function for exponential family;
    currently only binomial() with the logitistic linke is available in
    this version.}
  \item{n.models}{number of unique models to keep. If NULL, BAS will
    attempt to enumerate unless p > 35 or method="MCMC". For any of
    methods using MCMC algorithms that sample with replacement, sampling
    will stop when the number of iterations exceeds the min of
    'n.models' or 'MCMC.iterations' and on exit 'n.models' is updated to
    reflect the unique number of models that  have been sampled. }
   \item{betaprior}{Prior on coefficients for model coefficients (except
     intercept).  Options in clude CCH, robust, beta-prime, AIC, BIC. }
  \item{modelprior}{Family of prior distribution on the models.  Choices
    include \code{\link{uniform}}, \code{\link{Bernoulli}} or \code{\link{beta.binomial}.}}
  \item{initprobs}{vector of length p with the initial  inclusion
    probabilities used for sampling without replacement (the intercept
    will be included with probability one and does not need to be added here) or a character
    string giving the method used to construct the sampling probabilities
    if "Uniform" each predictor variable is equally likely to be
    sampled (equivalent to random sampling without replacement). If
    "eplogp", use the \code{\link{eplogprob}} function to aproximate the 
    Bayes factor using p-values to find initial marginal inclusion probabilitites and
    sample without replacement using these
    inclusion probabilaties, which may be updated using estimates of the
    marginal inclusion probabilites. "eplogp" assumes that MLEs from the
    full model exist; for problems where that is not the case or 'p' is
    large, initial sampling probabilities may be obtained using
    \code{\link{eplogprob.marg}} which fits a model to each predictor seaparately.
    For variables that should always be
    included set the corresponding initprobs to 1. To run a
    Markov Chain to provide initial estimates of marginal
    inclusion probabilities, use method="MCMC+BAS" below.}
  \item{method}{A character variable indicating which sampling method to
    use: method="BAS" uses Bayesian Adaptive Sampling (without
    replacement) using the sampling probabilities given in initprobs and
    updates using the marginal inclusion probabilities to direct the search/sample;
    method="MCMC" combines a random walk Metropolis Hastings (as in MC3 of
    Raftery et al 1997) with a random swap of a variable included with a
    variable that is currently excluded (see Clyde, Ghosh, and Littman
    (2010) for details);
    method="MCMC+BAS" runs an initial MCMC as above to calculate marginal
    inclusion probabilities and then samples without replacement as in
    BAS; method = "deterministic" runs an deterministic sampling using
    the initial probabilites (no updating); this is recommended for fast
    enumeration or if a model of independence is a good approximation to
    the joint posterior distribution of the model indicators.  For
    BAS, the sampling probabilities can be updated as more models are
    sampled. (see  'update' below).  We recommend "MCMC+BAS" or  "MCMC"
    for  high dimensional problems.} 
  \item{update}{number of iterations between potential updates of the
    sampling probabilities in the "BAS" method. If NULL do not update,
    otherwise the algorithm will update using the marginal inclusion
    probabilities as they change while sampling takes place.  For large
    model spaces, updating is recommended. If the model space will be
    enumerated,  leave at the default.}
  \item{bestmodel}{optional binary vector representing a model to
    initialize the sampling. If NULL sampling starts with the null
    model}
  \item{bestmarg}{optional value for the log marginal associated with
    the bestmodel}
  \item{prob.rw}{For any of the MCMC methods, probability of using the
    random-walk proposal; otherwise use a random "flip" move to propose
    a new model.}
  \item{MCMC.iterations}{Number of models to sample when using any of
    the MCMC options; should be greater than 'n.models'.}
  \item{control}{a list of parameters that control convergence in the
    fitting process.  See the documentation for
    \code{glm.control()}}
  \item{offset}{a priori known component to be included in the linear
predictor}
  \item{weights}{optional vector of weights to be used in the fitting
    process.  SHould be NULL or a numeric vector.}
  \item{laplace}{logical variable for whether to use a Laplace
    approximate for integration with respect to g to obtain the marginal
    likelihood.  If FALSE the Cephes library is used which may be
    inaccurate for large n or large values of the Wald Chisquared statistic.}
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
The deterinistic sampler provides a list of the top models in order of an
approximation of independence using the provided \code{initprobs}.  This
may be effective after running the other algorithms to identify high
probability models and works well if 
the correlations of variables are small to modest.  The priors on
coefficients are mixtures of g-priors that provide approximations to the
power prior.} 
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
  \item{Q}{the Q statistic for each model used in the 
    marginal likelihood approximation}
  \item{Y}{response}
  \item{X}{matrix of predictors}
}
\references{Li, Y. and Clyde, M. (2015) Mixtures of g-priors in
  Generalized Linear Models. \url{http://arxiv.org/abs/1503.06913}

 Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling
  for Variable Selection and Model Averaging. Journal of Computational
  Graphics and Statistics.  20:80-101 \cr \url{http://dx.doi.org/10.1198/jcgs.2010.09049}

  Raftery, A.E, Madigan, D. and Hoeting, J.A. (1997) Bayesian Model
  Averaging for Linear Regression Models. Journal of the American
  Statistical Association. 
}
\author{Merlise Clyde (\email{clyde@stat.duke.edu}), Quanli Wang and Yingbo Li}
\examples{
##---- Should be DIRECTLY executable !! ----
library(MASS)
data(Pima.tr)

out = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7, method="BAS",
 betaprior=CCH(a=1, b=532/2, s=0), family=binomial(),
 modelprior=beta.binomial(1,1), laplace=FALSE)

summary(out)
image(out)
}
\keyword{GLM}
\keyword{regression}
