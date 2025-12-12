# Fitted values for a BAS BMA objects

Calculate fitted values for a BAS BMA object

## Usage

``` r
# S3 method for class 'bas'
fitted(
  object,
  type = "link",
  estimator = "BMA",
  top = NULL,
  na.action = na.pass,
  ...
)
```

## Arguments

- object:

  An object of class 'bas' as created by
  [`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md)

- type:

  type equals "response" or "link" in the case of GLMs (default is
  'link')

- estimator:

  estimator type of fitted value to return. Default is to use BMA with
  all models. Options include  
  'HPM' the highest probability model  
  'BMA' Bayesian model averaging, using optionally only the 'top'
  models  
  'MPM' the median probability model of Barbieri and Berger. 'BPM' the
  model that is closest to BMA predictions under squared error loss

- top:

  optional argument specifying that the 'top' models will be used in
  constructing the BMA prediction, if NULL all models will be used. If
  top=1, then this is equivalent to 'HPM'

- na.action:

  function determining what should be done with missing values in
  newdata. The default is to predict NA.

- ...:

  optional arguments, not used currently

## Value

A vector of length n of fitted values.

## Details

Calculates fitted values at observed design matrix using either the
highest probability model, 'HPM', the posterior mean (under BMA) 'BMA',
the median probability model 'MPM' or the best predictive model 'BPM".
The median probability model is defined by including variable where the
marginal inclusion probability is greater than or equal to 1/2. For
type="BMA", the weighted average may be based on using a subset of the
highest probability models if an optional argument is given for top. By
default BMA uses all sampled models, which may take a while to compute
if the number of variables or number of models is large. The "BPM" is
found be computing the squared distance of the vector of fitted values
for a model and the fitted values under BMA and returns the model with
the smallest distance. In the presence of multicollinearity this may be
quite different from the MPM, with extreme collinearity may drop
relevant predictors.

## References

Barbieri, M. and Berger, J.O. (2004) Optimal predictive model selection.
Annals of Statistics. 32, 870-897.  
<https://projecteuclid.org/euclid.aos/1085408489&url=/UI/1.0/Summarize/euclid.aos/1085408489>

Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling
for Variable Selection and Model Averaging. Journal of Computational
Graphics and Statistics. 20:80-101  
[doi:10.1198/jcgs.2010.09049](https://doi.org/10.1198/jcgs.2010.09049)

## See also

[`predict.bas`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md)
[`predict.basglm`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md)

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/dev/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/coef.md),
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`diagnostics()`](http://merliseclyde.github.io/BAS/dev/reference/diagnostics.md),
[`force.heredity.bas()`](http://merliseclyde.github.io/BAS/dev/reference/force.heredity.bas.md),
[`image.bas()`](http://merliseclyde.github.io/BAS/dev/reference/image.bas.md),
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.confint.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/dev/reference/summary.md),
[`update.bas()`](http://merliseclyde.github.io/BAS/dev/reference/update.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

Other predict methods:
[`predict.bas()`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

## Author

Merlise Clyde <clyde@duke.edu>

## Examples

``` r
data(Hald)
hald.gprior =  bas.lm(Y~ ., data=Hald, prior="ZS-null", initprobs="Uniform")
plot(Hald$Y, fitted(hald.gprior, estimator="HPM"))

plot(Hald$Y, fitted(hald.gprior, estimator="BMA", top=3))

plot(Hald$Y, fitted(hald.gprior, estimator="MPM"))

plot(Hald$Y, fitted(hald.gprior, estimator="BPM"))

```
