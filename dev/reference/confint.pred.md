# Compute Credible (Bayesian Confidence) Intervals for a BAS predict object

Compute credible intervals for in-sample or out of sample prediction or
for the regression function

## Usage

``` r
# S3 method for class 'pred.bas'
confint(object, parm, level = 0.95, nsim = 10000, ...)
```

## Arguments

- object:

  an object created by
  [`predict.bas`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md)

- parm:

  character variable, "mean" or "pred". If missing parm='pred'.

- level:

  the nominal level of the (point-wise) credible interval

- nsim:

  number of Monte Carlo simulations for sampling methods with BMA

- ...:

  optional arguments to pass on to next function call; none at this
  time.

## Value

a matrix with lower and upper level \* 100 percent credible intervals
for either the mean of the regression function or predicted values.

## Details

This constructs approximate 95 percent Highest Posterior Density
intervals for 'pred.bas' objects. If the estimator is based on model
selection, the intervals use a Student t distribution using the estimate
of g. If the estimator is based on BMA, then nsim draws from the mixture
of Student t distributions are obtained with the HPD interval obtained
from the Monte Carlo draws.

## See also

[`predict.bas`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md)

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/dev/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/coef.md),
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.coef.md),
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

Other CI methods:
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.coef.md),
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.confint.md)

## Author

Merlise A Clyde

## Examples

``` r
data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
hald.pred = predict(hald.gprior, estimator="BPM", predict=FALSE, se.fit=TRUE)
confint(hald.pred, parm="mean")
#>            2.5%     97.5%      mean
#>  [1,]  72.85939  86.44363  79.65151
#>  [2,]  69.50215  79.45478  74.47846
#>  [3,] 101.93102 108.91264 105.42183
#>  [4,]  85.14928  94.51420  89.83174
#>  [5,]  89.98634 101.26964  95.62799
#>  [6,] 101.30948 107.88283 104.59616
#>  [7,]  97.81813 109.19556 103.50684
#>  [8,]  71.46153  82.55525  77.00839
#>  [9,]  87.73810  96.41333  92.07571
#> [10,] 106.51357 121.70394 114.10876
#> [11,]  77.32978  88.03488  82.68233
#> [12,] 106.79101 115.29472 111.04286
#> [13,] 105.64736 115.28746 110.46741
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"
confint(hald.pred)  #default
#>            2.5%     97.5%      pred
#>  [1,]  67.79843  91.50459  79.65151
#>  [2,]  63.56396  85.39296  74.47846
#>  [3,]  95.09960 115.74406 105.42183
#>  [4,]  79.04805 100.61543  89.83174
#>  [5,]  84.39452 106.86146  95.62799
#>  [6,]  94.34116 114.85115 104.59616
#>  [7,]  92.24966 114.76402 103.50684
#>  [8,]  65.82223  88.19456  77.00839
#>  [9,]  81.43722 102.71421  92.07571
#> [10,] 101.77792 126.43959 114.10876
#> [11,]  71.59123  93.77342  82.68233
#> [12,] 100.43905 121.64668 111.04286
#> [13,]  99.62326 121.31156 110.46741
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"
hald.pred = predict(hald.gprior, estimator="BMA", predict=FALSE, se.fit=TRUE)
confint(hald.pred)
#>            2.5%     97.5%      pred
#>  [1,]  67.54958  91.60154  79.68307
#>  [2,]  63.29096  85.83754  74.69127
#>  [3,]  94.98803 117.62923 105.63258
#>  [4,]  79.11098 101.38145  89.91648
#>  [5,]  84.03782 106.80151  95.67480
#>  [6,]  93.96497 114.94765 104.57616
#>  [7,]  92.05830 115.19454 103.47945
#>  [8,]  65.83565  88.43062  76.96808
#>  [9,]  81.26678 102.78682  92.22184
#> [10,] 100.99855 126.62569 113.84918
#> [11,]  71.77104  94.01569  82.59035
#> [12,]  99.86055 121.65736 110.87673
#> [13,]  99.06541 121.06843 110.34001
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"

```
