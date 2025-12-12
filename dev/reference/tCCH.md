# Generalized tCCH g-Prior Distribution for Coefficients in BMA Models

Creates an object representing the tCCH mixture of g-priors on
coefficients for BAS.

## Usage

``` r
tCCH(alpha = 1, beta = 2, s = 0, r = 3/2, v = 1, theta = 1)
```

## Arguments

- alpha:

  a scalar \> 0, recommended alpha=.5 (betaprime) or 1.

- beta:

  a scalar \> 0. The value is not updated by the data; beta should be a
  function of n for consistency under the null model.

- s:

  a scalar, recommended s=0 a priori

- r:

  r arbitrary; in the hyper-g-n prior sets r = (alpha + 2)

- v:

  0 \< v

- theta:

  theta \> 1

## Value

returns an object of class "prior", with the family and hyerparameters.

## Details

Creates a structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md).

## See also

[`CCH`](http://merliseclyde.github.io/BAS/dev/reference/CCH.md),
[`robust`](http://merliseclyde.github.io/BAS/dev/reference/robust.md),
[`hyper.g`](http://merliseclyde.github.io/BAS/dev/reference/hyper.g.md),
[`hyper.g.n`](http://merliseclyde.github.io/BAS/dev/reference/hyper.g.n.md)[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md)

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
[`robust()`](http://merliseclyde.github.io/BAS/dev/reference/robust.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/dev/reference/testBF.prior.md)

## Author

Merlise Clyde

## Examples

``` r
n <- 500
tCCH(alpha = 1, beta = 2, s = 0, r = 1.5, v = 1, theta = 1 / n)
#> $family
#> [1] "tCCH"
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
#> attr(,"class")
#> [1] "prior"
```
