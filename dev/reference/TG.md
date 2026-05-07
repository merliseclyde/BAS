# Generalized g-Prior Distribution for Coefficients in BMA Models

Creates an object representing the Truncated Gamma (tCCH) mixture of
g-priors on coefficients for BAS, where u = 1/(1+g) has a Gamma
distribution supported on (0, 1\].

## Usage

``` r
TG(alpha = 2)
```

## Arguments

- alpha:

  a scalar \> 0, recommended alpha=.5 (betaprime) or 1. alpha=2
  corresponds to the uniform prior on the shrinkage factor.

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md).

## See also

[`CCH`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md)
[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/dev/reference/EB.local.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/dev/reference/Jeffreys.md),
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

TG(alpha = 2)
#> $family
#> [1] "TG"
#> 
#> $class
#> [1] "TCCH"
#> 
#> $hyper.parameters
#> $hyper.parameters$alpha
#> [1] 2
#> 
#> $hyper.parameters$beta
#> [1] 2
#> 
#> $hyper.parameters$s
#> [1] 0
#> 
#> 
#> attr(,"class")
#> [1] "prior"
CCH(alpha = 2, beta = 100, s = 0)
#> $family
#> [1] "CCH"
#> 
#> $class
#> [1] "TCCH"
#> 
#> $hyper.parameters
#> $hyper.parameters$alpha
#> [1] 2
#> 
#> $hyper.parameters$beta
#> [1] 100
#> 
#> $hyper.parameters$s
#> [1] 0
#> 
#> 
#> attr(,"class")
#> [1] "prior"
```
