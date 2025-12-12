# Hyper-g-Prior Distribution for Coefficients in BMA Models

Creates an object representing the hyper-g mixture of g-priors on
coefficients for BAS.

## Usage

``` r
hyper.g(alpha = 3)
```

## Arguments

- alpha:

  a scalar \> 0. The hyper.g(alpha) is equivalent to CCH(alpha -2, 2,
  0). Liang et al recommended values in the range 2 \< alpha_h \<= 3

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
[`TG()`](http://merliseclyde.github.io/BAS/dev/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/dev/reference/beta.prime.md),
[`g.prior()`](http://merliseclyde.github.io/BAS/dev/reference/g.prior.md),
[`hyper.g.n()`](http://merliseclyde.github.io/BAS/dev/reference/hyper.g.n.md),
[`intrinsic()`](http://merliseclyde.github.io/BAS/dev/reference/intrinsic.md),
[`robust()`](http://merliseclyde.github.io/BAS/dev/reference/robust.md),
[`tCCH()`](http://merliseclyde.github.io/BAS/dev/reference/tCCH.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/dev/reference/testBF.prior.md)

## Author

Merlise Clyde

## Examples

``` r
hyper.g(alpha = 3)
#> $family
#> [1] "CCH"
#> 
#> $class
#> [1] "TCCH"
#> 
#> $hyper.parameters
#> $hyper.parameters$alpha
#> [1] 1
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
```
