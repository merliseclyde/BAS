# Post processing function to force constraints on interaction inclusion bas BMA objects

This function takes the output of a bas object and allows higher order
interactions to be included only if their parent lower order
interactions terms are in the model, by assigning zero prior
probability, and hence posterior probability, to models that do not
include their respective parents.

## Usage

``` r
force.heredity.bas(object, prior.prob = 0.5)
```

## Arguments

- object:

  a bas linear model or generalized linear model object

- prior.prob:

  prior probability that a term is included conditional on parents being
  included

## Value

a bas object with updated models, coefficients and summaries obtained
removing all models with zero prior and posterior probabilities.

## Note

Currently prior probabilities are computed using conditional Bernoulli
distributions, i.e. P(gamma_j = 1 \| Parents(gamma_j) = 1) = prior.prob.
This is not very efficient for models with a large number of levels.
Future updates will force this at the time of sampling.

## See also

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/dev/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/coef.md),
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`diagnostics()`](http://merliseclyde.github.io/BAS/dev/reference/diagnostics.md),
[`fitted.bas()`](http://merliseclyde.github.io/BAS/dev/reference/fitted.md),
[`image.bas()`](http://merliseclyde.github.io/BAS/dev/reference/image.bas.md),
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.confint.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/dev/reference/summary.md),
[`update.bas()`](http://merliseclyde.github.io/BAS/dev/reference/update.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

## Author

Merlise A Clyde

## Examples

``` r
data("chickwts")
bas.chk <- bas.lm(weight ~ feed, data = chickwts)
#  summary(bas.chk)  # 2^5 = 32 models
bas.chk.int <- force.heredity.bas(bas.chk)
#  summary(bas.chk.int)  # two models now


data(Hald)
bas.hald <- bas.lm(Y ~ .^2, data = Hald)
bas.hald.int <- force.heredity.bas(bas.hald)
image(bas.hald.int)


image(bas.hald.int)

# two-way interactions
data(ToothGrowth)
ToothGrowth$dose <- factor(ToothGrowth$dose)
levels(ToothGrowth$dose) <- c("Low", "Medium", "High")
TG.bas <- bas.lm(len ~ supp * dose, data = ToothGrowth, modelprior = uniform())
TG.bas.int <- force.heredity.bas(TG.bas)
image(TG.bas.int)
```
