# Plot Bayesian Confidence Intervals

Function takes the the output of functions that return credible
intervals from BAS objects, and creates a plot of the posterior mean
with segments representing the credible interval. of what the function
does. \~~

## Usage

``` r
# S3 method for class 'confint.bas'
plot(x, horizontal = FALSE, ...)
```

## Arguments

- x:

  the output from
  [`confint.coef.bas`](http://merliseclyde.github.io/BAS/reference/confint.coef.md)
  or
  [`confint.pred.bas`](http://merliseclyde.github.io/BAS/reference/confint.pred.md)
  containing credible intervals and estimates.

- horizontal:

  orientation of the plot

- ...:

  optional graphical arguments to pass on to plot

## Value

A plot of the credible intervals.

## Details

This function takes the HPD intervals or credible intervals created by
[`confint.coef.bas`](http://merliseclyde.github.io/BAS/reference/confint.coef.md)
or
[`confint.pred.bas`](http://merliseclyde.github.io/BAS/reference/confint.pred.md)
from BAS objects, and creates a plot of the posterior mean with segments
representing the credible interval. BAS tries to return HPD intervals,
and under model averaging these may not be symmetric. the description
above \~~

## See also

[`confint.coef.bas`](http://merliseclyde.github.io/BAS/reference/confint.coef.md),
[`confint.pred.bas`](http://merliseclyde.github.io/BAS/reference/confint.pred.md),
[`coef.bas`](http://merliseclyde.github.io/BAS/reference/coef.md),
[`predict.bas`](http://merliseclyde.github.io/BAS/reference/predict.bas.md),
`link{bas.lm}`

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/reference/coef.md),
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/reference/confint.coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/reference/confint.pred.md),
[`diagnostics()`](http://merliseclyde.github.io/BAS/reference/diagnostics.md),
[`fitted.bas()`](http://merliseclyde.github.io/BAS/reference/fitted.md),
[`force.heredity.bas()`](http://merliseclyde.github.io/BAS/reference/force.heredity.bas.md),
[`image.bas()`](http://merliseclyde.github.io/BAS/reference/image.bas.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/reference/summary.md),
[`update.bas()`](http://merliseclyde.github.io/BAS/reference/update.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/reference/variable.names.pred.bas.md)

Other CI methods:
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/reference/confint.coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/reference/confint.pred.md)

## Author

Merlise A Clyde

## Examples

``` r
data(Hald)
hald.ZS = bas.lm(Y ~ ., data=Hald, prior="ZS-null", modelprior=uniform())
hald.coef = confint(coef(hald.ZS), parm=2:5)
plot(hald.coef)

#> NULL
plot(hald.coef, horizontal=TRUE)

#> NULL
plot(confint(predict(hald.ZS, se.fit=TRUE), parm="mean"))

#> NULL
```
