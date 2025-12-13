# Families of G-Prior Distribution for Coefficients in BMA Models

Creates an object representing the g-prior distribution on coefficients
for BAS.

## Usage

``` r
g.prior(g)
```

## Arguments

- g:

  a scalar used in the covariance of Zellner's g-prior, Cov(beta) =
  sigma^2 g (X'X)^-1

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a structure used for BAS.

## See also

[`IC.prior`](http://merliseclyde.github.io/BAS/reference/IC.prior.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/reference/EB.local.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/reference/beta.prime.md),
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
g.prior(100)
#> $family
#> [1] "g.prior"
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
#> 
#> attr(,"class")
#> [1] "prior"
```
