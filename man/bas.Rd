\name{bas.lm}
\alias{bas}
\alias{bas.lm}
\title{Bayesian Adaptive Sampling Without Replacement for Variable
  Selection in Linear Models}
\description{Sample without replacement from a posterior distribution on models}
\usage{bas.lm(formula, data, weights=NULL, na.action="na.omit", 
 n.models=NULL,  prior="ZS-null", alpha=NULL,
 modelprior=beta.binomial(1,1),
 initprobs="Uniform", method="BAS", update=NULL,
 bestmodel = NULL, prob.local = 0.0, prob.rw=0.5,
 MCMC.iterations = NULL, lambda = NULL, delta = 0.025,thin=1)
}
\arguments{
  \item{formula}{linear model formula for the full model with all
    predictors, Y ~ X.  All code assumes that an intercept will be
    included in each model and that the X's will be centered.}
  \item{data}{data frame}
  \item{weights}{an optional vector of weights to be used in the fitting
    process. Should be NULL or a numeric vector. If non-NULL, Bayes
    estimates are obtained assuming that Y ~ N(Xb, sigma^2
    diag(1/weights)).}
    \item{na.action}{a function which indicates what should happen when the data contain NAs. The default is "na.omit".}
  \item{n.models}{number of models to sample either without replacement
    (method="BAS" or "MCMC+BAS")
    or with replacement (method="MCMC"). If NULL, BAS with
    method="BAS" will try to enumerate all 2^p models. If enumeration is
    not possible (memory or time) then a value should be supplied
    which controls the number of sampled models using 'n.models'.  
    With method="MCMC", sampling will stop once  at the min(n.models, MCMC.iterations) so MCMC.iterations  be larger than n.models in order to explore the model space. On exit for method=
    "MCMC"  this is the  number of unique models that have been
    sampled with counts stored in the output as "freq"" }
  \item{prior}{prior distribution for regression coefficients.  Choices include
    "AIC", "BIC", "g-prior", "ZS-null", "ZS-full", "hyper-g",
    "hyper-g-laplace", "hyper-g-n", "EB-local", and "EB-global" }
 \item{alpha}{optional hyperparameter in g-prior or hyper g-prior.  For
   Zellner's g-prior, alpha = g, for the Liang et al hyper-g or
   hyper-g-n method, recommended choice is alpha are between (2 < alpha < 4),
   with alpha = 3 recommended.} 
  \item{modelprior}{Family of prior distribution on the models.  Choices
    include \code{\link{uniform}}
    \code{\link{Bernoulli}} or \code{\link{beta.binomial}} with the
    default being a \code{beta.binomial(1,1)}.  }
  \item{initprobs}{Vector of length
    p or a character string for approach used to create the vector.
    This is used to order variables for sampling all methods for
    potentially more  efficient storage while sampling and provides the
    initial  inclusion  probabilities used for sampling without replacement with
    method="BAS".  Options for the charactier string giving the method
    are:  "Uniform" or "uniform" where each predictor variable is
    equally likely to be sampled (equivalent to random sampling without
    replacement);  "eplogp" uses the \code{\link{eplogprob}} function
    to aproximate the  Bayes factor from p-values from the full model to
    find initial marginal inclusion probabilitites;
    "marg-eplogp" uses\code{\link{eplogprob.marg}} function  to
    aproximate the  Bayes factor from p-values from the full model each
    simple linear regression.  To run a
    Markov Chain to provide initial estimates of marginal
    inclusion probabilities for "BAS", use method="MCMC+BAS" below.  While the
    initprobs are not used in sampling for method="MCMC", this
    determines the order of the variables in the lookup table and
    affects memory allocation in large problems where enumeration is not
    feasible.  For variables that should always be included set the
    corresponding initprobs to 1, i.e. the intercept should be included with
    probability one. 
   } 
  \item{method}{A character variable indicating which sampling method to
    use: method="BAS" uses Bayesian Adaptive Sampling (without
    replacement) using the sampling probabilities given in initprobs;
    method="MCMC" samples with replacement via a MCMC algorithm based
    that combines the birth/death random walk in Hoeting et al (1997) of
    MC3 with a random swap move to interchange a variable in the model
    with one currently excluded as in Clyde, Ghosh and Littman (2010);
    method="MCMC+BAS" runs an initial MCMC to calculate marginal
    inclusion probabilities and then samples without replacement as in
    BAS.  For both  BAS and AMCMC, the sampling probabilities can be updated as
    more  models are sampled. (see  update below).  We recommend
    "MCMC+BAS" or "MCMC" for high dimensional problems where enumeration
    is not feasible.}
  \item{update}{number of iterations between potential updates of the
    sampling probabilities for method "BAS". If NULL do not update,
    otherwise the  algorithm will update using the marginal inclusion
    probabilities as  they change while sampling takes place.  For large
    model spaces, updating is recommended. If the model space will be
    enumerated, leave at the default.}
  \item{bestmodel}{optional binary vector representing a model to
    initialize the sampling. If NULL sampling starts with the null
    model}
  \item{prob.local}{A future  option to allow sampling of models
    "near" the median probability model.  Not used at
    this time.}
  \item{prob.rw}{For any of the MCMC methods, probability of using the
    random-walk proposal; otherwise use a random "flip" move to propose
    a new model.}
  \item{MCMC.iterations}{Number of iterations for the MCMC sampler; the
    default is n.models*10 if not set by the user.} 
   \item{lambda}{Parameter in the AMCMC algorithm.}
   \item{delta}{truncation parameter to prevent sampling probabilities
     to degenerate to 0 or 1 prior to enumeration for sampling withour
     replacent.}
   \item{thin}{For "MCMC", thin the MCMC every "thin"
     iterations; default is no thinning.} 
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
coefficients include Zellner's g-prior, the Hyper-g prior (Liang et al
2008, the Zellner-Siow Cauchy prior, Empirical Bayes (local and gobal)
g-priors.  AIC and BIC are also included.  } \value{
  \code{bas} returns an object of class \code{bas}

An object of class \code{BAS} is a list containing at least the following components:

\item{postprob}{the posterior probabilities of the models selected}
\item{priorprobs}{the prior probabilities of the models selected}
  \item{namesx}{the names of the variables}
  \item{R2}{R2 values for the models}
  \item{logmarg}{values of the log of the marginal likelihood for the
    models}
 \item{n.vars}{total number of independent variables in the full model,
   including the intercept}
  \item{size}{the number of independent variables in each of the models,
    includes the intercept}
  \item{which}{a list of lists with one list per model with  variables
    that are included in the model} 
  \item{probne0}{the posterior probability that each variable is
    non-zero computed using the renormalized marginal likelihoods of
    sampled models.  This may be biased if the number of sampled models
    is much smaller than the total number of models. Unbiased estimates
    may be obtained using method "MCMC".}
  \item{mle}{list of lists with one list per model giving the MLE (OLS)
estimate of each (nonzero) coefficient for each model. NOTE: The intercept
is the mean of Y as each column of X has been centered by subtracting
its mean.}
  \item{mle.se}{list of lists with one list per model giving the MLE (OLS) standard error of each coefficient for each model}
  \item{prior}{the name of the prior  that created the BMA object}
  \item{alpha}{value of hyperparameter in prior used to create the BMA
    object.}
  \item{modelprior}{the prior distribution on models that created the BMA object}
  \item{Y}{response}
  \item{X}{matrix of predictors}
  \item{mean.x}{vector of means for each column of X (used in \code{\link{predict.bas}})}
  
 The function  \code{\link{summary.bas}}, is used to print a summary of
 the results. The function \code{\link{plot.bas}} is used to plot
 posterior distributions for the coefficients and
 \code{\link{image.bas}} provides an image of the distribution over models.
 Posterior summaries of coefficients can be extracted using
 \code{\link{coefficients.bas}}.  Fitted values and predictions can be
 obtained using the S3 functions  \code{\link{fitted.bas}} and \code{\link{predict.bas}}.
 BAS objects may be updated to use a different prior (without rerunning
 the sampler) using the function  \code{\link{update.bas}}.
}

\references{
  Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling
  for Variable Selection and Model Averaging. Journal of Computational
  Graphics and Statistics.  20:80-101 \cr \url{http://dx.doi.org/10.1198/jcgs.2010.09049}
  
  Clyde, M. and George, E. I. (2004) Model Uncertainty. Statist. Sci.,
  19, 81-94. \cr \url{http://dx.doi.org/10.1214/088342304000000035}
  
  Clyde, M. (1999)
  Bayesian Model Averaging and Model Search Strategies (with
  discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid,
  J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages
  157-185.

  Hoeting, J. A., Madigan, D., Raftery, A. E. and Volinsky, C. T. (1999)
  Bayesian model averaging: a tutorial (with discussion). Statist. Sci.,
  14, 382-401. \cr \url{http://www.stat.washington.edu/www/research/online/hoeting1999.pdf}
  
  Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger,
  J.O. (2008) Mixtures of  g-priors for Bayesian Variable
  Selection. Journal of the American Statistical Association.
   103:410-423.  \cr \url{http://dx.doi.org/10.1198/016214507000001337}

  
  Zellner, A. (1986) On assessing prior distributions and Bayesian
  regression analysis with g-prior distributions. In Bayesian Inference
  and Decision Techniques: Essays in Honor of Bruno de Finetti,
  pp. 233-243. North-Holland/Elsevier.
    
  Zellner, A. and Siow, A. (1980) Posterior odds ratios for selected
  regression hypotheses. In Bayesian Statistics: Proceedings of the First
  International Meeting held in Valencia (Spain), pp. 585-603.
}


\author{Merlise Clyde (\email{clyde@stat.duke.edu}) and Michael Littman}
  
\seealso{ \code{\link{summary.bas}},
          \code{\link{coefficients.bas}},
          \code{\link{print.bas}},
	  \code{\link{predict.bas}},
	  \code{\link{fitted.bas}}
          \code{\link{plot.bas}},
	  \code{\link{image.bas}},
	  \code{\link{eplogprob}},
	  \code{\link{update.bas}}
	  
}
	
\examples{
demo(BAS.hald)
\dontrun{demo(BAS.USCrime) }
}

\keyword{regression}


