# Empirical Bayes Prior Distribution for Coefficients in BMA Model

Creates an object representing the EB prior for BAS GLM.

## Usage

``` r
EB.local()
```

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/reference/bas.glm.md).

## See also

[`CCH`](http://merliseclyde.github.io/BAS/reference/CCH.md) and
[`bas.glm`](http://merliseclyde.github.io/BAS/reference/bas.glm.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/reference/CCH.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/reference/beta.prime.md),
[`g.prior()`](http://merliseclyde.github.io/BAS/reference/g.prior.md),
[`hyper.g()`](http://merliseclyde.github.io/BAS/reference/hyper.g.md),
[`hyper.g.n()`](http://merliseclyde.github.io/BAS/reference/hyper.g.n.md),
[`intrinsic()`](http://merliseclyde.github.io/BAS/reference/intrinsic.md),
[`robust()`](http://merliseclyde.github.io/BAS/reference/robust.md),
[`tCCH()`](http://merliseclyde.github.io/BAS/reference/tCCH.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/reference/testBF.prior.md)

## Author

Merlise Clyde

## Examples

``` r
EB.local()
#> $family
#> [1] "EB-local"
#> 
#> $class
#> [1] "EB"
#> 
#> $hyper.parameters
#> $hyper.parameters$local
#> [1] TRUE
#> 
#> 
#> attr(,"class")
#> [1] "prior"
```
