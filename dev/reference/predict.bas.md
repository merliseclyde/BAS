# Prediction Method for an object of class BAS

Predictions under model averaging or other estimators from a BMA object
of class inheriting from 'bas'.

## Usage

``` r
# S3 method for class 'bas'
predict(
  object,
  newdata,
  se.fit = FALSE,
  type = c("response", "link"),
  top = NULL,
  estimator = "BMA",
  na.action = na.pass,
  ...
)
```

## Arguments

- object:

  An object of class BAS, created by `bas`

- newdata:

  dataframe for predictions. If missing, then use the dataframe used for
  fitting for obtaining fitted and predicted values.

- se.fit:

  indicator for whether to compute se of fitted and predicted values

- type:

  Type of predictions required. The defaults "reponse" is on the scale
  of the response is the only option currently for linear models (for
  Gaussian models this is equivalent to type="link").

- top:

  a scalar integer M. If supplied, subset the top M models, based on
  posterior probabilities for model predictions and BMA.

- estimator:

  estimator used for predictions. Currently supported options include:  
  'HPM' the highest probability model  
  'BMA' Bayesian model averaging, using optionally only the 'top'
  models  
  'MPM' the median probability model of Barbieri and Berger.  
  'BPM' the model that is closest to BMA predictions under squared error
  loss. BMA may be computed using only the 'top' models if supplied

- na.action:

  function determining what should be done with missing values in
  newdata. The default is to predict NA.

- ...:

  optional extra arguments

## Value

a list of

- fit:

  fitted values based on the selected estimator

- Ybma:

  predictions using BMA, the same as fit for non-BMA methods for
  compatibility; will be deprecated

- Ypred:

  matrix of predictions under each model for BMA

- se.fit:

  se of fitted values; in the case of BMA this will be a matrix

- se.pred:

  se for predicted values; in the case of BMA this will be a matrix

- se.bma.fit:

  vector of posterior sd under BMA for posterior mean of the regression
  function. This will be NULL if estimator is not 'BMA'

- se.bma.pred:

  vector of posterior sd under BMA for posterior predictive values. this
  will be NULL if estimator is not 'BMA'

- best:

  index of top models included

- bestmodels:

  subset of bestmodels used for fitting or prediction

- best.vars:

  names of variables in the top model; NULL if estimator='BMA'

- df:

  scalar or vector of degrees of freedom for models

- estimator:

  estimator upon which 'fit' is based.

## Details

Use BMA and/or model selection to form predictions using the top highest
probability models.

## See also

[`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`fitted.bas`](http://merliseclyde.github.io/BAS/dev/reference/fitted.md),
[`confint.pred.bas`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`variable.names.pred.bas`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

Other predict methods:
[`fitted.bas()`](http://merliseclyde.github.io/BAS/dev/reference/fitted.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

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
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/dev/reference/summary.md),
[`update.bas()`](http://merliseclyde.github.io/BAS/dev/reference/update.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

## Author

Merlise Clyde

## Examples

``` r

data("Hald")
hald.gprior =  bas.lm(Y ~ ., data=Hald, alpha=13, prior="g-prior")

predict(hald.gprior, newdata=Hald, estimator="BPM", se.fit=TRUE)
#> $fit
#>  [1]  79.65151  74.47846 105.42183  89.83174  95.62799 104.59616 103.50684
#>  [8]  77.00839  92.07571 114.10876  82.68233 111.04286 110.46741
#> attr(,"model")
#> [1] 0 1 2 4
#> attr(,"best")
#> [1] 14
#> attr(,"estimator")
#> [1] "BPM"
#> 
#> $Ybma
#>            [,1]
#>  [1,]  79.68307
#>  [2,]  74.69127
#>  [3,] 105.63258
#>  [4,]  89.91648
#>  [5,]  95.67480
#>  [6,] 104.57616
#>  [7,] 103.47945
#>  [8,]  76.96808
#>  [9,]  92.22184
#> [10,] 113.84918
#> [11,]  82.59035
#> [12,] 110.87673
#> [13,] 110.34001
#> 
#> $Ypred
#>            [,1]     [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,]  81.17036 74.83464 105.07248  89.69881  97.15898 104.45753 103.38927
#>  [2,]  77.70296 74.24113 105.85537  90.46267  93.09565 104.71517 103.13993
#>  [3,]  79.70437 74.40553 105.21752  89.76253  95.63309 104.57088 103.52541
#>  [4,]  79.65151 74.47846 105.42183  89.83174  95.62799 104.59616 103.50684
#>  [5,]  79.84321 74.31409 104.90632  89.65651  95.70301 104.52849 103.54760
#>  [6,]  79.26248 74.75042 106.28314  90.16732  95.37830 104.70862 103.39893
#>  [7,]  78.80123 75.69457 108.22151  90.62061  95.54467 104.84276 103.50032
#>  [8,]  81.66556 77.02098 106.35099  88.18423  99.83232 103.89115 105.74346
#>  [9,]  85.78065 79.39070 104.28069  87.30339 103.43704 102.66524 106.03984
#> [10,]  74.86000 80.34349 102.27744  83.77067  93.36677 100.90656 111.87354
#> [11,]  79.18965 81.38793 101.17241  82.85345  98.24138 100.43966 112.16380
#> [12,]  76.29823 80.55874 101.93127  83.25765  95.26054 100.79397 112.20198
#> [13,]  94.62219 84.21059 101.56325 101.56325  94.62219 101.56325  87.68112
#> [14,]  95.42308 95.42308  95.42308  95.42308  95.42308  95.42308  95.42308
#> [15,]  91.78308 83.03167 101.29055 101.29055  91.78308 101.74970  88.24455
#> [16,] 102.15048 91.65573  99.81831  99.81831 102.15048  98.65223  89.32357
#>           [,8]      [,9]     [,10]    [,11]     [,12]     [,13]
#>  [1,] 76.06454  91.57174 113.17222 81.59906 111.22195 111.08841
#>  [2,] 78.80193  92.68123 115.80581 84.50293 110.41616 109.07906
#>  [3,] 77.08557  91.98604 114.17593 82.78145 111.11959 110.53210
#>  [4,] 77.00839  92.07571 114.10876 82.68233 111.04286 110.46741
#>  [5,] 77.15919  91.83513 114.23530 82.88128 111.23840 110.65147
#>  [6,] 76.86019  92.49134 113.99208 82.44826 110.67744 110.08148
#>  [7,] 76.13446  93.59932 112.64184 81.53107 109.86900 109.49863
#>  [8,] 74.60469  93.86382 106.77052 80.21897 110.61958 111.73373
#>  [9,] 74.19437  93.55892 101.91430 79.36984 110.13525 112.42979
#> [10,] 85.82697 100.90656  98.16482 92.68133 107.76092 107.76092
#> [11,] 82.85345  99.70690  94.57759 89.44827 108.50000 109.96552
#> [12,] 84.53056 100.50527  96.78718 91.37187 108.21267 108.79006
#> [13,] 84.21059  85.94586 118.91590 84.21059 101.56325  99.82798
#> [14,] 95.42308  95.42308  95.42308 95.42308  95.42308  95.42308
#> [15,] 86.24572  86.55641 120.92687 86.70486 101.74970  99.14326
#> [16,] 83.49316  88.15749 104.48264 82.32707  98.65223  99.81831
#> 
#> $postprobs
#>  [1] 2.432256e-01 1.684081e-01 1.312165e-01 1.224293e-01 1.220358e-01
#>  [6] 1.145513e-01 6.888252e-02 2.709377e-02 1.481347e-03 2.891971e-04
#> [11] 2.559516e-04 5.597869e-05 4.790052e-05 1.177702e-05 1.000762e-05
#> [16] 5.009504e-06
#> 
#> $se.fit
#>        1        2        3        4        5        6        7        8 
#> 3.117350 2.283957 1.602160 2.149087 2.589321 1.508471 2.610923 2.545817 
#>        9       10       11       12       13 
#> 1.990817 3.485929 2.456636 1.951456 2.212238 
#> 
#> $se.pred
#>        1        2        3        4        5        6        7        8 
#> 5.440156 5.009380 4.737547 4.949344 5.155775 4.706689 5.166658 5.134065 
#>        9       10       11       12       13 
#> 4.882702 5.659428 5.090431 4.866787 4.977090 
#> 
#> $se.bma.fit
#> NULL
#> 
#> $se.bma.pred
#> NULL
#> 
#> $df
#> [1] 12
#> 
#> $best
#> [1] 14
#> 
#> $bestmodel
#> [1] 0 1 2 4
#> 
#> $best.vars
#> [1] "Intercept" "X1"        "X2"        "X4"       
#> 
#> $estimator
#> [1] "BPM"
#> 
#> attr(,"class")
#> [1] "pred.bas"
# same as fitted
fitted(hald.gprior,estimator="BPM")
#>  [1]  79.65151  74.47846 105.42183  89.83174  95.62799 104.59616 103.50684
#>  [8]  77.00839  92.07571 114.10876  82.68233 111.04286 110.46741
# default is BMA and estimation of mean vector
hald.bma = predict(hald.gprior, top=5, se.fit=TRUE)
confint(hald.bma)
#>            2.5%     97.5%      pred
#>  [1,]  67.81695  91.35950  79.74246
#>  [2,]  63.52803  85.44744  74.50010
#>  [3,]  94.61178 116.49325 105.29268
#>  [4,]  78.61249 100.47169  89.88693
#>  [5,]  85.24122 107.10752  95.57177
#>  [6,]  94.44972 115.27340 104.56409
#>  [7,]  92.25425 115.22190 103.40145
#>  [8,]  65.57941  88.30339  77.13668
#>  [9,]  81.32509 102.88603  91.99731
#> [10,] 101.39902 126.65998 114.21325
#> [11,]  70.99091  93.43341  82.78446
#> [12,] 100.54377 121.80017 111.00723
#> [13,]  99.13533 121.62011 110.40160
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"

hald.bpm = predict(hald.gprior, newdata=Hald[1,],
                    se.fit=TRUE,
                    estimator="BPM")
confint(hald.bpm)
#>          2.5%    97.5%     pred
#> [1,] 67.74454 91.66421 79.70437
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"
# extract variables
variable.names(hald.bpm)
#> [1] "Intercept" "X1"        "X2"        "X3"        "X4"       

hald.hpm = predict(hald.gprior, newdata=Hald[1,],
                    se.fit=TRUE,
                    estimator="HPM")
confint(hald.hpm)
#>          2.5%    97.5%     pred
#> [1,] 70.15171 92.18902 81.17036
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"
variable.names(hald.hpm)
#> [1] "Intercept" "X1"        "X2"       

hald.mpm = predict(hald.gprior, newdata=Hald[1,],
                    se.fit=TRUE,
                    estimator="MPM")
confint(hald.mpm)
#>          2.5%    97.5%     pred
#> [1,] 67.79843 91.50459 79.65151
#> attr(,"Probability")
#> [1] 0.95
#> attr(,"class")
#> [1] "confint.bas"
variable.names(hald.mpm)
#> [1] "Intercept" "X1"        "X2"        "X4"       
```
