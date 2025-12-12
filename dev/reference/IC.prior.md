# Information Criterion Families of Prior Distribution for Coefficients in BMA Models

Creates an object representing the prior distribution on coefficients
for BAS.

## Usage

``` r
IC.prior(penalty)
```

## Arguments

- penalty:

  a scalar used in the penalized loglikelihood of the form
  penalty\*dimension

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

The log marginal likelihood is approximated as -2\*(deviance +
penalty\*dimension). Allows alternatives to AIC (penalty = 2) and BIC
(penalty = log(n)). For BIC, the argument may be missing, in which case
the sample size is determined from the call to \`bas.glm\` and used to
determine the penalty.

## See also

[`g.prior`](http://merliseclyde.github.io/BAS/dev/reference/g.prior.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/dev/reference/EB.local.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/dev/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/dev/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/dev/reference/beta.prime.md),
[`g.prior()`](http://merliseclyde.github.io/BAS/dev/reference/g.prior.md),
[`hyper.g()`](http://merliseclyde.github.io/BAS/dev/reference/hyper.g.md),
[`hyper.g.n()`](http://merliseclyde.github.io/BAS/dev/reference/hyper.g.n.md),
[`intrinsic()`](http://merliseclyde.github.io/BAS/dev/reference/intrinsic.md),
[`robust()`](http://merliseclyde.github.io/BAS/dev/reference/robust.md),
[`tCCH()`](http://merliseclyde.github.io/BAS/dev/reference/tCCH.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/dev/reference/testBF.prior.md)

## Author

Merlise Clyde

## Examples

``` r
IC.prior(2)
#> $family
#> [1] "IC"
#> 
#> $class
#> [1] "IC"
#> 
#> $hyper
#> [1] 2
#> 
#> $hyper.parameters
#> $hyper.parameters$penalty
#> [1] 2
#> 
#> 
#> attr(,"class")
#> [1] "prior"
aic.prior()
#> $family
#> [1] "AIC"
#> 
#> $class
#> [1] "IC"
#> 
#> $hyper.parameters
#> $hyper.parameters$penalty
#> [1] 2
#> 
#> 
#> $hyper
#> [1] 2
#> 
#> attr(,"class")
#> [1] "prior"
bic.prior(100)
#> $family
#> [1] "BIC"
#> 
#> $class
#> [1] "IC"
#> 
#> $hyper.parameters
#> $hyper.parameters$penalty
#> [1] 4.60517
#> 
#> $hyper.parameters$n
#> [1] 100
#> 
#> 
#> $hyper
#> [1] 4.60517
#> 
#> attr(,"class")
#> [1] "prior"
```
