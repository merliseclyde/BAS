# BAS 1.3.1   July 26, 2016

## New features

* new function to plot credible intervals created by `confint.pred.bas` or `confint.coef.bas`.   See the help files for an example or the vignette.


# BAS 1.3.0   July 15, 2016

## New Features

A vignette has been added at long last!  This illustrates several of the new features in `BAS` such as
  
* new functions for computing credible intervals for fitted and predicted values  `confint.pred.bas()`
* new function for adding credible intervals for coefficients  `confint.coef.bas()` 
* added posterior standard deviations for fitted values and predicted values in  `predict.bas()`
  


## Deprication 
* deprecated use of `type` to specify estimator in fitted.bas	and replaced with `estimator` so that `predict()` and `fitted()` are compatible with other S3 methods. 
* updated funtions to be of class `bas` to avoid NAMESPACE conficts with  other libraries
	

# BAS 1.2.2 June 29, 2016

## New Features
* added option to find "Best Predictive Model" or "BPM" for `fitted.bas` or `predict.bas`
* added local Empirical Bayes prior and fixed g-prior for `bas.glm`
* added `diagnostic()` function for checking convergence of `bas` objects created with `method = "MCMC"`"
* added truncated power prior as in Yang, Wainwright & Jordan (2016)

##	Minor Changes

* bug fix in `plot.bas` that appears with Sweave
* bug fix in `coef.bma` when there is just one predictor
	
	

# BAS 1.2.1  April 16, 2016
* bug fix for method="MCMC" with truncated prior distributions
	where MH ratio was incorrect allowing models with 0 probability to
	be sampled.
* fixed error in Zellner-Siow prior (ZS-null) when n=p+1 or
	saturated model  where log marginal likelihood should be 0 

# BAS 1.2.0  April 11, 2016
* removed unsafe code where Rbestmarg (input) was being
	overwritten in .Call which would end up in corruption of the
	constant pool of the byte-code  (Thanks to Tomas Kalibera for
	catching this!)
* fixed issue with dimensions for use with Simple Linear Regression

# BAS 1.1.0   March 31, 2016

## New Features
* added truncated Beta-Binomial prior and truncated Poisson (works
	only with MCMC currently)
* improved code for finding fitted values under the Median 
* deprecated method = "AMCMC" and issue warning message

## Minor Changes
* Changed S3 method for plot and image to use class `bas` rather than
	 `bma` to avoid name conflicts with other packages

# BAS 1.09
	- added weights for linear models
	- switched LINPACK calls in bayesreg to LAPACK finally should be
	faster
	- fixed bug in intercept calculation for glms
	- fixed inclusion probabilities to be a vector in the global EB
	methods for linear models
# BAS 1.08
	- added intrinsic prior for GLMs
	- fixed problems for linear models for p > n and R2 not correct
# BAS 1.07
	- added phi1 function from Gordy (1998)  confluent hypergeometric
	function of two variables  also known as one of the Horn
	hypergeometric functions or Humbert's phi1
	- added Jeffrey's prior on g
	- added the general tCCH prior and special cases of the hyper-g/n.
	- TODO check shrinkage functions for all	
# BAS 1.06
	- new improved Laplace approximation for hypergeometric1F1
	- added class basglm for predict
	- predict function now handles glm output
	- added dataframe option for newdata in predict.bas and predict.basglm
	- renamed coefficients in output to be 'mle' in bas.lm to be consistent across
	lm and glm versions so that predict methods can handle both
	cases.  (This may lead to errors in other external code that
	expects object$ols or object$coefficients)
	- fixed bug with initprobs that did not include an intercept for bas.lm
	
# BAS 1.05
	- added thinning option for MCMC method for bas.lm
	- returned posterior expected shrinkage for bas.glm
	- added option for initprobs = "marg-eplogp" for using marginal
	SLR models to create starting probabilities or order variables
	especially for p > n case
	- added standalone function for hypergeometric1F1 using Cephes
	library and a Laplace aproximation
	-Added class "BAS" so that predict and fitted functions (S3
	methods) are not masked by functions in the BVS package: to do
	modify the rest of the S3 methods.
	
# BAS 1.04

	- added bas.glm for model averaging/section using mixture of g-priors for
	GLMs.  Currently limited to Logistic Regression
	- added Poisson family for glm.fit

# BAS 1.0	
	- cleaned up  MCMC method code
	
# BAS 0.93

	- removed internal print statements in bayesglm.c
	- Bug fixes in AMCMC algorithm

# BAS 0.92

	- fixed glm-fit.R  so that hyperparameter for BIC is numeric

# BAS 0.91

	- added new AMCMC algorithm

# BAS 0.91

	- bug fix in bayes.glm

# BAS 0.90

	- added C routines for fitting glms

# BAS 0.85

	- fixed problem with duplicate models if n.models was > 2^(p-1) by
   restricting n.models

	- save original X as part of object so that fitted.bma gives the
   correct fitted values (broken in version 0.80)
 
# BAS 0.80

	- Added `hypergeometric2F1` function that is callable by R
	- centered X's in bas.lm so that the intercept has the correct
  shrinkage
	- changed `predict.bma` to center newdata using the mean(X)
	- Added new Adaptive MCMC option (method = "AMCMC")  (this is not
  stable at this point)

# BAS 0.7

	-Allowed pruning of model tree to eliminate rejected models
 
# BAS 0.6

	- Added MCMC option to create starting values for BAS (`method = "MCMC+BAS"`)

# BAS 0.5

	-Cleaned up all .Call routines so that all objects are duplicated or
 allocated within code

# BAS 0.45

	- fixed ch2inv that prevented building on Windows in bayes glm_fit

# BAS 0.4

	- fixed fortran calls to use F77_NAME macro 
	- changed  allocation of objects for .Call to prevent some objects from being overwritten.  
  
# BAS 0.3

	- fixed EB.global function to include prior probabilities on models
	- fixed update function 

# BAS 0.2

	- fixed predict.bma to allow newdata to be a matrix or vector with the
  column of ones for the intercept optionally included.
	- fixed help file for predict 
	- added modelprior argument to bas.lm so that users may now use the
	beta-binomial prior distribution on model size in addition to the
	default uniform distribution
	- added functions uniform(), beta-binomial() and Bernoulli() to create
	model prior objects
	- added a vector of user specified initial probabilities as an option for
	argument initprobs in bas.lm and removed the separate argument user.prob

