# Package index

## model fitting

Functions associated with inference and selection for Bayesian Linear
and Generalized Linear Models

- [`bas.glm()`](http://merliseclyde.github.io/BAS/reference/bas.glm.md)
  : Bayesian Adaptive Sampling Without Replacement for Variable
  Selection in Generalized Linear Models
- [`bas.lm()`](http://merliseclyde.github.io/BAS/reference/bas.lm.md) :
  Bayesian Adaptive Sampling for Bayesian Model Averaging and Variable
  Selection in Linear Models
- [`bayesglm.fit()`](http://merliseclyde.github.io/BAS/reference/bayesglm.fit.md)
  : Fitting Generalized Linear Models and Bayesian marginal likelihood
  evaluation
- [`coef(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/coef.md)
  [`print(`*`<coef.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/coef.md)
  : Coefficients of a Bayesian Model Average object
- [`confint(`*`<coef.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/confint.coef.md)
  : Compute Credible Intervals for BAS regression coefficients from BAS
  objects
- [`confint(`*`<pred.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/confint.pred.md)
  : Compute Credible (Bayesian Confidence) Intervals for a BAS predict
  object
- [`cv.summary.bas()`](http://merliseclyde.github.io/BAS/reference/cv.summary.bas.md)
  : Summaries for Out of Sample Prediction
- [`diagnostics()`](http://merliseclyde.github.io/BAS/reference/diagnostics.md)
  : BAS MCMC diagnostic plot
- [`eplogprob()`](http://merliseclyde.github.io/BAS/reference/eplogprob.md)
  : eplogprob - Compute approximate marginal inclusion probabilities
  from pvalues
- [`eplogprob.marg()`](http://merliseclyde.github.io/BAS/reference/eplogprob.marg.md)
  : eplogprob.marg - Compute approximate marginal inclusion
  probabilities from pvalues
- [`fitted(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/fitted.md)
  : Fitted values for a BAS BMA objects
- [`force.heredity.bas()`](http://merliseclyde.github.io/BAS/reference/force.heredity.bas.md)
  : Post processing function to force constraints on interaction
  inclusion bas BMA objects
- [`image(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/image.bas.md)
  : Images of models used in Bayesian model averaging
- [`plot(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/plot.md)
  : Plot Diagnostics for an BAS Object
- [`plot(`*`<coef.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/plot.coef.md)
  : Plots the posterior distributions of coefficients derived from
  Bayesian model averaging
- [`plot(`*`<confint.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/plot.confint.md)
  : Plot Bayesian Confidence Intervals
- [`predict(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/predict.bas.md)
  : Prediction Method for an object of class BAS
- [`predict(`*`<basglm>`*`)`](http://merliseclyde.github.io/BAS/reference/predict.basglm.md)
  : Prediction Method for an Object of Class basglm
- [`print(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/print.bas.md)
  : Print a Summary of Bayesian Model Averaging objects from BAS
- [`summary(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/summary.md)
  : Summaries of Bayesian Model Averaging objects from BAS
- [`update(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/update.md)
  : Update BAS object using a new prior
- [`variable.names(`*`<pred.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/variable.names.pred.bas.md)
  : Extract the variable names for a model from a BAS prediction object

## coefficient priors

- [`beta.prime()`](http://merliseclyde.github.io/BAS/reference/beta.prime.md)
  : Beta-Prime Prior Distribution for Coefficients in BMA Model
- [`CCH()`](http://merliseclyde.github.io/BAS/reference/CCH.md) :
  Generalized g-Prior Distribution for Coefficients in BMA Models
- [`EB.global()`](http://merliseclyde.github.io/BAS/reference/EB.global.md)
  : Find the global Empirical Bayes estimates for BMA
- [`EB.local()`](http://merliseclyde.github.io/BAS/reference/EB.local.md)
  : Empirical Bayes Prior Distribution for Coefficients in BMA Model
- [`g.prior()`](http://merliseclyde.github.io/BAS/reference/g.prior.md)
  : Families of G-Prior Distribution for Coefficients in BMA Models
- [`hyper.g()`](http://merliseclyde.github.io/BAS/reference/hyper.g.md)
  : Hyper-g-Prior Distribution for Coefficients in BMA Models
- [`hyper.g.n()`](http://merliseclyde.github.io/BAS/reference/hyper.g.n.md)
  : Generalized hyper-g/n Prior Distribution for g for mixtures of
  g-priors on Coefficients in BMA Models
- [`IC.prior()`](http://merliseclyde.github.io/BAS/reference/IC.prior.md)
  : Information Criterion Families of Prior Distribution for
  Coefficients in BMA Models
- [`intrinsic()`](http://merliseclyde.github.io/BAS/reference/intrinsic.md)
  : Intrinsic Prior Distribution for Coefficients in BMA Models
- [`Jeffreys()`](http://merliseclyde.github.io/BAS/reference/Jeffreys.md)
  : Jeffreys Prior Distribution for \$g\$ for Mixtures of g-Priors for
  Coefficients in BMA Models
- [`robust()`](http://merliseclyde.github.io/BAS/reference/robust.md) :
  Robust-Prior Distribution for Coefficients in BMA Model
- [`tCCH()`](http://merliseclyde.github.io/BAS/reference/tCCH.md) :
  Generalized tCCH g-Prior Distribution for Coefficients in BMA Models
- [`TG()`](http://merliseclyde.github.io/BAS/reference/TG.md) :
  Generalized g-Prior Distribution for Coefficients in BMA Models
- [`testBF.prior()`](http://merliseclyde.github.io/BAS/reference/testBF.prior.md)
  : Test based Bayes Factors for BMA Models

## model priors

- [`Bernoulli()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md)
  : Independent Bernoulli Prior Distribution for Models
- [`Bernoulli.heredity()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.heredity.md)
  : Independent Bernoulli prior on models that with constraints for
  model hierarchy induced by interactions
- [`beta.binomial()`](http://merliseclyde.github.io/BAS/reference/beta.binomial.md)
  : Beta-Binomial Prior Distribution for Models
- [`tr.beta.binomial()`](http://merliseclyde.github.io/BAS/reference/tr.beta.binomial.md)
  : Truncated Beta-Binomial Prior Distribution for Models
- [`tr.poisson()`](http://merliseclyde.github.io/BAS/reference/tr.poisson.md)
  : Truncated Poisson Prior Distribution for Models
- [`tr.power.prior()`](http://merliseclyde.github.io/BAS/reference/tr.power.prior.md)
  : Truncated Power Prior Distribution for Models
- [`uniform()`](http://merliseclyde.github.io/BAS/reference/uniform.md)
  : Uniform Prior Distribution for Models

## special functions

- [`hypergeometric1F1()`](http://merliseclyde.github.io/BAS/reference/hypergeometric1F1.md)
  : Confluent hypergeometric1F1 function
- [`hypergeometric2F1()`](http://merliseclyde.github.io/BAS/reference/hypergeometric2F1.md)
  : Gaussian hypergeometric2F1 function
- [`phi1()`](http://merliseclyde.github.io/BAS/reference/phi1.md) :
  Compound Confluent hypergeometric function of two variables

## utility functions

- [`Bayes.outlier()`](http://merliseclyde.github.io/BAS/reference/Bayes.outlier.md)
  : Bayesian Outlier Detection
- [`list2matrix.bas()`](http://merliseclyde.github.io/BAS/reference/list2matrix.md)
  : Coerce a BAS list object into a matrix.
- [`list2matrix.which()`](http://merliseclyde.github.io/BAS/reference/list2matrix.which.md)
  : Coerce a BAS list object into a matrix.
- [`which.matrix()`](http://merliseclyde.github.io/BAS/reference/which.matrix.md)
  : Coerce a BAS list object of models into a matrix.

## data sets

- [`Hald`](http://merliseclyde.github.io/BAS/reference/Hald.md)
  [`hald`](http://merliseclyde.github.io/BAS/reference/Hald.md) : Hald
  Data
- [`bodyfat`](http://merliseclyde.github.io/BAS/reference/bodyfat.md)
  [`Bodyfat`](http://merliseclyde.github.io/BAS/reference/bodyfat.md) :
  Bodyfat Data
- [`protein`](http://merliseclyde.github.io/BAS/reference/protein.md) :
  Protein Activity Data

## all functions

- [`BAS`](http://merliseclyde.github.io/BAS/reference/BAS.md) : BAS:
  Bayesian Model Averaging using Bayesian Adaptive Sampling
- [`Bayes.outlier()`](http://merliseclyde.github.io/BAS/reference/Bayes.outlier.md)
  : Bayesian Outlier Detection
- [`Bernoulli()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md)
  : Independent Bernoulli Prior Distribution for Models
- [`Bernoulli.heredity()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.heredity.md)
  : Independent Bernoulli prior on models that with constraints for
  model hierarchy induced by interactions
- [`CCH()`](http://merliseclyde.github.io/BAS/reference/CCH.md) :
  Generalized g-Prior Distribution for Coefficients in BMA Models
- [`EB.global()`](http://merliseclyde.github.io/BAS/reference/EB.global.md)
  : Find the global Empirical Bayes estimates for BMA
- [`EB.local()`](http://merliseclyde.github.io/BAS/reference/EB.local.md)
  : Empirical Bayes Prior Distribution for Coefficients in BMA Model
- [`Hald`](http://merliseclyde.github.io/BAS/reference/Hald.md)
  [`hald`](http://merliseclyde.github.io/BAS/reference/Hald.md) : Hald
  Data
- [`IC.prior()`](http://merliseclyde.github.io/BAS/reference/IC.prior.md)
  : Information Criterion Families of Prior Distribution for
  Coefficients in BMA Models
- [`Jeffreys()`](http://merliseclyde.github.io/BAS/reference/Jeffreys.md)
  : Jeffreys Prior Distribution for \$g\$ for Mixtures of g-Priors for
  Coefficients in BMA Models
- [`TG()`](http://merliseclyde.github.io/BAS/reference/TG.md) :
  Generalized g-Prior Distribution for Coefficients in BMA Models
- [`bas.glm()`](http://merliseclyde.github.io/BAS/reference/bas.glm.md)
  : Bayesian Adaptive Sampling Without Replacement for Variable
  Selection in Generalized Linear Models
- [`bas.lm()`](http://merliseclyde.github.io/BAS/reference/bas.lm.md) :
  Bayesian Adaptive Sampling for Bayesian Model Averaging and Variable
  Selection in Linear Models
- [`bayesglm.fit()`](http://merliseclyde.github.io/BAS/reference/bayesglm.fit.md)
  : Fitting Generalized Linear Models and Bayesian marginal likelihood
  evaluation
- [`beta.binomial()`](http://merliseclyde.github.io/BAS/reference/beta.binomial.md)
  : Beta-Binomial Prior Distribution for Models
- [`beta.prime()`](http://merliseclyde.github.io/BAS/reference/beta.prime.md)
  : Beta-Prime Prior Distribution for Coefficients in BMA Model
- [`bodyfat`](http://merliseclyde.github.io/BAS/reference/bodyfat.md)
  [`Bodyfat`](http://merliseclyde.github.io/BAS/reference/bodyfat.md) :
  Bodyfat Data
- [`climate`](http://merliseclyde.github.io/BAS/reference/climate.md) :
  Climate Data
- [`coef(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/coef.md)
  [`print(`*`<coef.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/coef.md)
  : Coefficients of a Bayesian Model Average object
- [`confint(`*`<coef.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/confint.coef.md)
  : Compute Credible Intervals for BAS regression coefficients from BAS
  objects
- [`confint(`*`<pred.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/confint.pred.md)
  : Compute Credible (Bayesian Confidence) Intervals for a BAS predict
  object
- [`cv.summary.bas()`](http://merliseclyde.github.io/BAS/reference/cv.summary.bas.md)
  : Summaries for Out of Sample Prediction
- [`diagnostics()`](http://merliseclyde.github.io/BAS/reference/diagnostics.md)
  : BAS MCMC diagnostic plot
- [`eplogprob()`](http://merliseclyde.github.io/BAS/reference/eplogprob.md)
  : eplogprob - Compute approximate marginal inclusion probabilities
  from pvalues
- [`eplogprob.marg()`](http://merliseclyde.github.io/BAS/reference/eplogprob.marg.md)
  : eplogprob.marg - Compute approximate marginal inclusion
  probabilities from pvalues
- [`fitted(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/fitted.md)
  : Fitted values for a BAS BMA objects
- [`force.heredity.bas()`](http://merliseclyde.github.io/BAS/reference/force.heredity.bas.md)
  : Post processing function to force constraints on interaction
  inclusion bas BMA objects
- [`g.prior()`](http://merliseclyde.github.io/BAS/reference/g.prior.md)
  : Families of G-Prior Distribution for Coefficients in BMA Models
- [`hyper.g()`](http://merliseclyde.github.io/BAS/reference/hyper.g.md)
  : Hyper-g-Prior Distribution for Coefficients in BMA Models
- [`hyper.g.n()`](http://merliseclyde.github.io/BAS/reference/hyper.g.n.md)
  : Generalized hyper-g/n Prior Distribution for g for mixtures of
  g-priors on Coefficients in BMA Models
- [`hypergeometric1F1()`](http://merliseclyde.github.io/BAS/reference/hypergeometric1F1.md)
  : Confluent hypergeometric1F1 function
- [`hypergeometric2F1()`](http://merliseclyde.github.io/BAS/reference/hypergeometric2F1.md)
  : Gaussian hypergeometric2F1 function
- [`image(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/image.bas.md)
  : Images of models used in Bayesian model averaging
- [`intrinsic()`](http://merliseclyde.github.io/BAS/reference/intrinsic.md)
  : Intrinsic Prior Distribution for Coefficients in BMA Models
- [`list2matrix.bas()`](http://merliseclyde.github.io/BAS/reference/list2matrix.md)
  : Coerce a BAS list object into a matrix.
- [`list2matrix.which()`](http://merliseclyde.github.io/BAS/reference/list2matrix.which.md)
  : Coerce a BAS list object into a matrix.
- [`phi1()`](http://merliseclyde.github.io/BAS/reference/phi1.md) :
  Compound Confluent hypergeometric function of two variables
- [`plot(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/plot.md)
  : Plot Diagnostics for an BAS Object
- [`plot(`*`<coef.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/plot.coef.md)
  : Plots the posterior distributions of coefficients derived from
  Bayesian model averaging
- [`plot(`*`<confint.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/plot.confint.md)
  : Plot Bayesian Confidence Intervals
- [`predict(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/predict.bas.md)
  : Prediction Method for an object of class BAS
- [`predict(`*`<basglm>`*`)`](http://merliseclyde.github.io/BAS/reference/predict.basglm.md)
  : Prediction Method for an Object of Class basglm
- [`print(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/print.bas.md)
  : Print a Summary of Bayesian Model Averaging objects from BAS
- [`protein`](http://merliseclyde.github.io/BAS/reference/protein.md) :
  Protein Activity Data
- [`robust()`](http://merliseclyde.github.io/BAS/reference/robust.md) :
  Robust-Prior Distribution for Coefficients in BMA Model
- [`summary(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/summary.md)
  : Summaries of Bayesian Model Averaging objects from BAS
- [`tCCH()`](http://merliseclyde.github.io/BAS/reference/tCCH.md) :
  Generalized tCCH g-Prior Distribution for Coefficients in BMA Models
- [`testBF.prior()`](http://merliseclyde.github.io/BAS/reference/testBF.prior.md)
  : Test based Bayes Factors for BMA Models
- [`tr.beta.binomial()`](http://merliseclyde.github.io/BAS/reference/tr.beta.binomial.md)
  : Truncated Beta-Binomial Prior Distribution for Models
- [`tr.poisson()`](http://merliseclyde.github.io/BAS/reference/tr.poisson.md)
  : Truncated Poisson Prior Distribution for Models
- [`tr.power.prior()`](http://merliseclyde.github.io/BAS/reference/tr.power.prior.md)
  : Truncated Power Prior Distribution for Models
- [`trCCH()`](http://merliseclyde.github.io/BAS/reference/trCCH.md) :
  Truncated Compound Confluent Hypergeometric function
- [`uniform()`](http://merliseclyde.github.io/BAS/reference/uniform.md)
  : Uniform Prior Distribution for Models
- [`update(`*`<bas>`*`)`](http://merliseclyde.github.io/BAS/reference/update.md)
  : Update BAS object using a new prior
- [`variable.names(`*`<pred.bas>`*`)`](http://merliseclyde.github.io/BAS/reference/variable.names.pred.bas.md)
  : Extract the variable names for a model from a BAS prediction object
- [`which.matrix()`](http://merliseclyde.github.io/BAS/reference/which.matrix.md)
  : Coerce a BAS list object of models into a matrix.
