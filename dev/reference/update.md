# Update BAS object using a new prior

Update a BMA object using a new prior distribution on the coefficients.

## Usage

``` r
# S3 method for class 'bas'
update(object, newprior, alpha = NULL, ...)
```

## Arguments

- object:

  BMA object to update

- newprior:

  Update posterior model probabilities, probne0, shrinkage, logmarg,
  etc, using prior based on newprior. See
  [`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md) for
  available methods

- alpha:

  optional new value of hyperparameter in prior for method

- ...:

  optional arguments

## Value

A new object of class BMA

## Details

Recomputes the marginal likelihoods for the new methods for models
already sampled in current object.

## References

Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling
for Variable Selection and Model Averaging. Journal of Computational
Graphics and Statistics. 20:80-101  
[doi:10.1198/jcgs.2010.09049](https://doi.org/10.1198/jcgs.2010.09049)

## See also

[`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md) for
available methods and choices of alpha

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/dev/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/coef.md),
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`diagnostics()`](http://merliseclyde.github.io/BAS/dev/reference/diagnostics.md),
[`fitted.bas()`](http://merliseclyde.github.io/BAS/dev/reference/fitted.md),
[`force.heredity.bas()`](http://merliseclyde.github.io/BAS/dev/reference/force.heredity.bas.md),
[`image.bas()`](http://merliseclyde.github.io/BAS/dev/reference/image.bas.md),
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.confint.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/dev/reference/summary.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

## Author

Merlise Clyde <clyde@stat.duke.edu>

## Examples

``` r

# \donttest{
library(MASS)
data(UScrime)
UScrime[,-2] <- log(UScrime[,-2])
crime.bic <-  bas.lm(y ~ ., data=UScrime, n.models=2^10, prior="BIC",initprobs= "eplogp")
crime.ebg <- update(crime.bic, newprior="EB-global")
crime.zs <- update(crime.bic, newprior="ZS-null")
# }
```
