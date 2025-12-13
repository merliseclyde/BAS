# Summaries of Bayesian Model Averaging objects from BAS

`summary` and `print` methods for Bayesian model averaging objects
created by `bas` Bayesian Adaptive Sampling

## Usage

``` r
# S3 method for class 'bas'
summary(object, n.models = 5, ...)
```

## Arguments

- object:

  object of class 'bas'

- n.models:

  optional number specifying the number of best models to display in
  summary

- ...:

  other parameters to be passed to `summary.default`

## Details

The print methods display a view similar to `print.lm` . The summary
methods display a view specific to Bayesian model averaging giving the
top 5 highest probability models represented by their inclusion
indicators. Summaries of the models include the Bayes Factor (BF) of
each model to the model with the largest marginal likelihood, the
posterior probability of the models, R2, dim (which includes the
intercept) and the log of the marginal likelihood.

## See also

[`coef.bas`](http://merliseclyde.github.io/BAS/reference/coef.md)

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
[`update.bas()`](http://merliseclyde.github.io/BAS/reference/update.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/reference/variable.names.pred.bas.md)

## Author

Merlise Clyde <clyde@duke.edu>

## Examples

``` r
data(UScrime, package = "MASS")
UScrime[, -2] <- log(UScrime[, -2])
crime.bic <- bas.lm(y ~ ., data = UScrime, n.models = 2^15, prior = "BIC", initprobs = "eplogp")
print(crime.bic)
#> 
#> Call:
#> bas.lm(formula = y ~ ., data = UScrime, n.models = 2^15, prior = "BIC", 
#>     initprobs = "eplogp")
#> 
#> 
#>  Marginal Posterior Inclusion Probabilities: 
#> Intercept          M         So         Ed        Po1        Po2         LF  
#>    1.0000     0.9335     0.3277     0.9910     0.7247     0.4602     0.2935  
#>       M.F        Pop         NW         U1         U2        GDP       Ineq  
#>    0.3298     0.4963     0.8346     0.3481     0.7752     0.5254     0.9992  
#>      Prob       Time  
#>    0.9541     0.5433  
summary(crime.bic)
#>           P(B != 0 | Y)   model 1       model 2     model 3     model 4
#> Intercept     1.0000000   1.00000  1.000000e+00   1.0000000   1.0000000
#> M             0.9335117   1.00000  1.000000e+00   1.0000000   1.0000000
#> So            0.3276563   0.00000  1.000000e+00   0.0000000   0.0000000
#> Ed            0.9910219   1.00000  1.000000e+00   1.0000000   1.0000000
#> Po1           0.7246635   1.00000  1.000000e+00   1.0000000   1.0000000
#> Po2           0.4602481   0.00000  1.000000e+00   0.0000000   0.0000000
#> LF            0.2935326   0.00000  1.000000e+00   0.0000000   0.0000000
#> M.F           0.3298168   0.00000  1.000000e+00   0.0000000   0.0000000
#> Pop           0.4962869   0.00000  1.000000e+00   0.0000000   0.0000000
#> NW            0.8346412   1.00000  1.000000e+00   1.0000000   1.0000000
#> U1            0.3481266   0.00000  1.000000e+00   0.0000000   0.0000000
#> U2            0.7752102   1.00000  1.000000e+00   1.0000000   1.0000000
#> GDP           0.5253694   0.00000  1.000000e+00   0.0000000   1.0000000
#> Ineq          0.9992058   1.00000  1.000000e+00   1.0000000   1.0000000
#> Prob          0.9541470   1.00000  1.000000e+00   1.0000000   1.0000000
#> Time          0.5432686   1.00000  1.000000e+00   0.0000000   1.0000000
#> BF                   NA   1.00000  1.267935e-04   0.7609295   0.5431578
#> PostProbs            NA   0.01910  1.560000e-02   0.0145000   0.0133000
#> R2                   NA   0.84200  8.695000e-01   0.8265000   0.8506000
#> dim                  NA   9.00000  1.600000e+01   8.0000000  10.0000000
#> logmarg              NA -22.15855 -3.113150e+01 -22.4317627 -22.7689035
#>               model 5
#> Intercept   1.0000000
#> M           1.0000000
#> So          0.0000000
#> Ed          1.0000000
#> Po1         1.0000000
#> Po2         0.0000000
#> LF          0.0000000
#> M.F         0.0000000
#> Pop         1.0000000
#> NW          1.0000000
#> U1          0.0000000
#> U2          1.0000000
#> GDP         0.0000000
#> Ineq        1.0000000
#> Prob        1.0000000
#> Time        0.0000000
#> BF          0.5203179
#> PostProbs   0.0099000
#> R2          0.8375000
#> dim         9.0000000
#> logmarg   -22.8118635
```
