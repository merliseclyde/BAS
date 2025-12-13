# Test based Bayes Factors for BMA Models

Creates an object representing the prior distribution on coefficients
for BAS that corresponds to the test-based Bayes Factors.

## Usage

``` r
testBF.prior(g)
```

## Arguments

- g:

  a scalar used in the covariance of Zellner's g-prior, Cov(beta) =
  sigma^2 g (X'X)^-

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a prior object structure used for BAS in \`bas.glm\`.

## See also

[`g.prior`](http://merliseclyde.github.io/BAS/reference/g.prior.md),
[`bas.glm`](http://merliseclyde.github.io/BAS/reference/bas.glm.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/reference/EB.local.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/reference/beta.prime.md),
[`g.prior()`](http://merliseclyde.github.io/BAS/reference/g.prior.md),
[`hyper.g()`](http://merliseclyde.github.io/BAS/reference/hyper.g.md),
[`hyper.g.n()`](http://merliseclyde.github.io/BAS/reference/hyper.g.n.md),
[`intrinsic()`](http://merliseclyde.github.io/BAS/reference/intrinsic.md),
[`robust()`](http://merliseclyde.github.io/BAS/reference/robust.md),
[`tCCH()`](http://merliseclyde.github.io/BAS/reference/tCCH.md)

## Author

Merlise Clyde

## Examples

``` r
testBF.prior(100)
#> $family
#> [1] "testBF.prior"
#> 
#> $g
#> [1] 100
#> 
#> $class
#> [1] "g-prior"
#> 
#> $hyper
#> [1] 100
#> 
#> $hyper.parameters
#> $hyper.parameters$g
#> [1] 100
#> 
#> $hyper.parameters$loglik_null
#> NULL
#> 
#> 
#> attr(,"class")
#> [1] "prior"
library(MASS)
data(Pima.tr)

# use g = n
bas.glm(type ~ .,
  data = Pima.tr, family = binomial(),
  betaprior = testBF.prior(nrow(Pima.tr)),
  modelprior = uniform(), method = "BAS"
)
#> 
#> Call:
#> bas.glm(formula = type ~ ., family = binomial(), data = Pima.tr, 
#>     betaprior = testBF.prior(nrow(Pima.tr)), modelprior = uniform(), 
#>     method = "BAS")
#> 
#> 
#>  Marginal Posterior Inclusion Probabilities: 
#> Intercept      npreg        glu         bp       skin        bmi        ped  
#>    1.0000     0.4252     1.0000     0.0706     0.1264     0.6139     0.8075  
#>       age  
#>    0.6705  
```
