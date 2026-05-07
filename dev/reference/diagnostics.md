# BAS MCMC diagnostic plot

Function to help assess convergence of MCMC sampling for bas objects.

## Usage

``` r
diagnostics(obj, type = c("pip", "model"), ...)
```

## Arguments

- obj:

  an object created by bas.lm or bas.glm

- type:

  type of diagnostic plot. If "pip" the marginal inclusion probabilities
  are used, while if "model", plot posterior model probabilities

- ...:

  additional graphics parameters to be passed to plot

## Value

a plot with of the marginal inclusion probabilities (pip) estimated by
MCMC and renormalized marginal likelihoods times prior probabilities or
model probabilities.

## Details

BAS calculates posterior model probabilities in two ways when
method="MCMC". The first is using the relative Monte Carlo frequencies
of sampled models. The second is to renormalize the marginal likelihood
times prior probabilities over the sampled models. If the Markov chain
has converged, these two quantities should be the same and fall on a 1-1
line. If not, running longer may be required. If the chain has not
converged, the Monte Carlo frequencies may have less bias, although may
exhibit more variability on repeated runs.

## See also

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/dev/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/coef.md),
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`fitted.bas()`](http://merliseclyde.github.io/BAS/dev/reference/fitted.md),
[`force.heredity.bas()`](http://merliseclyde.github.io/BAS/dev/reference/force.heredity.bas.md),
[`image.bas()`](http://merliseclyde.github.io/BAS/dev/reference/image.bas.md),
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.confint.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/dev/reference/summary.md),
[`update.bas()`](http://merliseclyde.github.io/BAS/dev/reference/update.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

## Author

Merlise Clyde (<clyde@duke.edu>)

## Examples

``` r

library(MASS)
data(UScrime)
UScrime[, -2] <- log(UScrime[, -2])
crime.ZS <- bas.lm(y ~ .,
  data = UScrime,
  prior = "ZS-null",
  modelprior = uniform(),
  method = "MCMC",
  MCMC.iter = 1000
) # short run for the example
diagnostics(crime.ZS)

```
