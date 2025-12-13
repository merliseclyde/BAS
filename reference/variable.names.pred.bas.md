# Extract the variable names for a model from a BAS prediction object

S3 method for class 'pred.bas'. Simple utility function to extract the
variable names. Used to print names for the selected models using
estimators for 'HPM', 'MPM' or 'BPM". for the selected model created by
`predict` for BAS objects.

## Usage

``` r
# S3 method for class 'pred.bas'
variable.names(object, ...)
```

## Arguments

- object:

  a BAS object created by `predict` from a BAS \`bas.lm\` or \`bas.glm\`
  object

- ...:

  other arguments to pass on

## Value

a character vector with the names of the variables included in the
selected model; in the case of 'BMA' this will be all variables

## See also

[`predict.bas`](http://merliseclyde.github.io/BAS/reference/predict.bas.md)

Other predict methods:
[`fitted.bas()`](http://merliseclyde.github.io/BAS/reference/fitted.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/reference/predict.basglm.md)

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
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/reference/plot.confint.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/reference/summary.md),
[`update.bas()`](http://merliseclyde.github.io/BAS/reference/update.md)

## Examples

``` r
data(Hald)
hald.gprior =  bas.lm(Y~ ., data=Hald, prior="ZS-null", modelprior=uniform())
hald.bpm = predict(hald.gprior, newdata=Hald[1,],
                   se.fit=TRUE,
                   estimator="BPM")
variable.names(hald.bpm)
#> [1] "Intercept" "X2"       
```
