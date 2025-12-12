# Beta-Prime Prior Distribution for Coefficients in BMA Model

Creates an object representing the Beta-Prime prior that is mixture of
g-priors on coefficients for BAS.

## Usage

``` r
beta.prime(n = NULL)
```

## Arguments

- n:

  the sample size; if NULL, the value derived from the data in the call
  to \`bas.glm\` will be used.

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md).

## See also

[`CCH`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/dev/reference/EB.local.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/dev/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/dev/reference/TG.md),
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
beta.prime(n = 100)
#> $family
#> [1] "betaprime"
#> 
#> $class
#> [1] "TCCH"
#> 
#> $hyper.parameters
#> $hyper.parameters$n
#> [1] 100
#> 
#> $hyper.parameters$alpha
#> [1] 0.5
#> 
#> 
#> attr(,"class")
#> [1] "prior"
```
