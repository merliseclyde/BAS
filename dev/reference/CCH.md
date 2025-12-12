# Generalized g-Prior Distribution for Coefficients in BMA Models

Creates an object representing the CCH mixture of g-priors on
coefficients for BAS .

## Usage

``` r
CCH(alpha, beta, s = 0)
```

## Arguments

- alpha:

  a scalar \> 0, recommended alpha=.5 (betaprime) or 1 for CCH. The
  hyper.g(alpha) is equivalent to CCH(alpha -2, 2, 0). Liang et al
  recommended values in the range 2 \< alpha_h \<= 4

- beta:

  a scalar \> 0. The value is not updated by the data; beta should be a
  function of n for consistency under the null model. The hyper-g
  corresponds to b = 2

- s:

  a scalar, recommended s=0

## Value

returns an object of class "prior", with the family and hyperparameters.

## Details

Creates a structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md).

## See also

[`IC.prior`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md),
[`bic.prior`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md),
[`bas.glm`](http://merliseclyde.github.io/BAS/dev/reference/bas.glm.md)

Other beta priors:
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
[`tCCH()`](http://merliseclyde.github.io/BAS/dev/reference/tCCH.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/dev/reference/testBF.prior.md)

## Author

Merlise A Clyde

## Examples

``` r
CCH(alpha = .5, beta = 100, s = 0)
#> $family
#> [1] "CCH"
#> 
#> $class
#> [1] "TCCH"
#> 
#> $hyper.parameters
#> $hyper.parameters$alpha
#> [1] 0.5
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
