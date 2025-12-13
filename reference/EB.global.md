# Find the global Empirical Bayes estimates for BMA

Finds the global Empirical Bayes estimates of g in Zellner's g-prior and
model probabilities

## Usage

``` r
EB.global(object, tol = 0.1, g.0 = NULL, max.iterations = 100)
```

## Arguments

- object:

  A 'bas' object created by
  [`bas`](http://merliseclyde.github.io/BAS/reference/bas.lm.md)

- tol:

  tolerance for estimating g

- g.0:

  initial value for g

- max.iterations:

  Maximum number of iterations for the EM algorithm

## Value

An object of class 'bas' using Zellner's g prior with an estimate of g
based on all models

## Details

Uses the EM algorithm in Liang et al to estimate the type II MLE of g in
Zellner's g prior

## References

Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2008)
Mixtures of g-priors for Bayesian Variable Selection. Journal of the
American Statistical Association. 103:410-423.  
[doi:10.1198/016214507000001337](https://doi.org/10.1198/016214507000001337)

## See also

[`bas`](http://merliseclyde.github.io/BAS/reference/bas.lm.md),
[`update`](http://merliseclyde.github.io/BAS/reference/update.md)

## Author

Merlise Clyde <clyde@stat.duke.edu>

## Examples

``` r
library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
# EB local uses a different g within each model
crime.EBL =  bas.lm(y ~ ., data=UScrime, n.models=2^15,
                    prior="EB-local", initprobs= "eplogp")
# use a common (global) estimate of g
crime.EBG = EB.global(crime.EBL)

```
