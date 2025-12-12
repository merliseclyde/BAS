# Robust-Prior Distribution for Coefficients in BMA Model

Creates an object representing the robust prior of Bayarri et al (2012)
that is mixture of g-priors on coefficients for BAS.

## Usage

``` r
robust(n = NULL)
```

## Arguments

- n:

  the sample size.

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a prior structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md).

## See also

[`CCH`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md)
and[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/dev/reference/EB.local.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/dev/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/dev/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/dev/reference/beta.prime.md),
[`g.prior()`](http://merliseclyde.github.io/BAS/dev/reference/g.prior.md),
[`hyper.g()`](http://merliseclyde.github.io/BAS/dev/reference/hyper.g.md),
[`hyper.g.n()`](http://merliseclyde.github.io/BAS/dev/reference/hyper.g.n.md),
[`intrinsic()`](http://merliseclyde.github.io/BAS/dev/reference/intrinsic.md),
[`tCCH()`](http://merliseclyde.github.io/BAS/dev/reference/tCCH.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/dev/reference/testBF.prior.md)

## Author

Merlise Clyde

## Examples

``` r
robust(100)
#> $family
#> [1] "robust"
#> 
#> $class
#> [1] "TCCH"
#> 
#> $hyper.parameters
#> $hyper.parameters$n
#> [1] 100
#> 
#> 
#> attr(,"class")
#> [1] "prior"
```
