# Compute Credible Intervals for BAS regression coefficients from BAS objects

Uses Monte Carlo simulations using posterior means and standard
deviations of coefficients to generate draws from the posterior
distributions and returns highest posterior density (HPD) credible
intervals. If the number of models equals one, then use the t
distribution to find intervals. These currently condition on the
estimate of \$g\$. than the description above \~~

## Usage

``` r
# S3 method for class 'coef.bas'
confint(object, parm, level = 0.95, nsim = 10000, ...)
```

## Arguments

- object:

  a coef.bas object

- parm:

  a specification of which parameters are to be given credible
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  the probability coverage required

- nsim:

  number of Monte Carlo draws from the posterior distribution. Used when
  number of models is greater than 1.

- ...:

  other arguments to passed; none currently

## Value

A matrix (or vector) with columns giving lower and upper HPD credible
limits for each parameter. These will be labeled as 1-level)/2 and 1 -
(1-level)/2 in percent (by default 2.5 and 97.5).

## Note

For mixture of g-priors these are approximate. This uses Monte Carlo
sampling so results may be subject to Monte Carlo variation and larger
values of nsim may be needed to reduce variability.

## See also

Other CI methods:
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.confint.md)

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/dev/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`diagnostics()`](http://merliseclyde.github.io/BAS/dev/reference/diagnostics.md),
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

Merlise A Clyde

## Examples

``` r

data("Hald")
hald_gprior <-  bas.lm(Y~ ., data=Hald, alpha=13,
                            prior="g-prior")
coef_hald <- coef(hald_gprior)
confint(coef_hald)
#>                2.5%      97.5%       beta
#> Intercept 93.826357 96.9017427 95.4230769
#> X1         0.000000  1.8732079  1.2150202
#> X2        -1.059551  0.9839873  0.2756235
#> X3        -1.511040  0.5923170 -0.1270575
#> X4        -1.660959  0.3338077 -0.3268710
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"
confint(coef_hald, approx=FALSE, nsim=5000)
#>                2.5%      97.5%       beta
#> Intercept 93.842978 96.8887411 95.4230769
#> X1         0.000000  1.8809801  1.2150202
#> X2        -1.103061  0.8990128  0.2756235
#> X3        -1.522862  0.5544236 -0.1270575
#> X4        -1.693533  0.2582418 -0.3268710
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"
# extract just the coefficient of X4
confint(coef_hald, parm="X4")
#>         2.5%     97.5%      beta
#> X4 -1.670591 0.2483164 -0.326871
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"

```
