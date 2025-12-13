# Generalized hyper-g/n Prior Distribution for g for mixtures of g-priors on Coefficients in BMA Models

Creates an object representing the hyper-g/n mixture of g-priors on
coefficients for BAS. This is a special case of the tCCH prior

## Usage

``` r
hyper.g.n(alpha = 3, n = NULL)
```

## Arguments

- alpha:

  a scalar \> 0, recommended 2 \< alpha \<= 3

- n:

  The sample size; if NULL, the value derived from the data in the call
  to \`bas.glm\` will be used.

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/reference/bas.glm.md).
This is a special case of the
[`tCCH`](http://merliseclyde.github.io/BAS/reference/tCCH.md), where
`hyper.g.n(alpha=3, n)` is equivalent to
` tCCH(alpha=1, beta=2, s=0, r=1.5, v = 1, theta=1/n) `

## See also

[`tCCH`](http://merliseclyde.github.io/BAS/reference/tCCH.md),
[`robust`](http://merliseclyde.github.io/BAS/reference/robust.md),
[`hyper.g`](http://merliseclyde.github.io/BAS/reference/hyper.g.md),
[`CCH`](http://merliseclyde.github.io/BAS/reference/CCH.md)[`bas.glm`](http://merliseclyde.github.io/BAS/reference/bas.glm.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/reference/EB.local.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/reference/beta.prime.md),
[`g.prior()`](http://merliseclyde.github.io/BAS/reference/g.prior.md),
[`hyper.g()`](http://merliseclyde.github.io/BAS/reference/hyper.g.md),
[`intrinsic()`](http://merliseclyde.github.io/BAS/reference/intrinsic.md),
[`robust()`](http://merliseclyde.github.io/BAS/reference/robust.md),
[`tCCH()`](http://merliseclyde.github.io/BAS/reference/tCCH.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/reference/testBF.prior.md)

## Author

Merlise Clyde

## Examples

``` r
n <- 500
hyper.g.n(alpha = 3, n = n)
#> $family
#> [1] "hyper-g/n"
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
#> $hyper.parameters$r
#> [1] 1.5
#> 
#> $hyper.parameters$v
#> [1] 1
#> 
#> $hyper.parameters$theta
#> [1] 0.002
#> 
#> 
#> $n
#> [1] 500
#> 
#> attr(,"class")
#> [1] "prior"
```
