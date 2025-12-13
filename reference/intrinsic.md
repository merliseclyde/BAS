# Intrinsic Prior Distribution for Coefficients in BMA Models

Creates an object representing the intrinsic prior on g, a special case
of the tCCH mixture of g-priors on coefficients for BAS.

## Usage

``` r
intrinsic(n = NULL)
```

## Arguments

- n:

  the sample size; if NULL, the value derived from the data in the call
  to \`bas.glm\` will be used.

## Value

returns an object of class "prior", with the family "intrinsic" of class
"TCCH" and hyperparameters alpha = 1, beta = 1, s = 0, r = 1, n = n for
the tCCH prior where theta in the tCCH prior is determined by the model
size and sample size.

## Details

Creates a structure used for
[`bas.glm`](http://merliseclyde.github.io/BAS/reference/bas.glm.md).

## References

Womack, A., Novelo,L.L., Casella, G. (2014). "Inference From Intrinsic
Bayes' Procedures Under Model Selection and Uncertainty". Journal of the
American Statistical Association. 109:1040-1053.
[doi:10.1080/01621459.2014.880348](https://doi.org/10.1080/01621459.2014.880348)

## See also

[`tCCH`](http://merliseclyde.github.io/BAS/reference/tCCH.md),
[`robust`](http://merliseclyde.github.io/BAS/reference/robust.md),
[`hyper.g`](http://merliseclyde.github.io/BAS/reference/hyper.g.md),
[`hyper.g.n`](http://merliseclyde.github.io/BAS/reference/hyper.g.n.md)[`bas.glm`](http://merliseclyde.github.io/BAS/reference/bas.glm.md)

Other beta priors:
[`CCH()`](http://merliseclyde.github.io/BAS/reference/CCH.md),
[`EB.local()`](http://merliseclyde.github.io/BAS/reference/EB.local.md),
[`IC.prior()`](http://merliseclyde.github.io/BAS/reference/IC.prior.md),
[`Jeffreys()`](http://merliseclyde.github.io/BAS/reference/Jeffreys.md),
[`TG()`](http://merliseclyde.github.io/BAS/reference/TG.md),
[`beta.prime()`](http://merliseclyde.github.io/BAS/reference/beta.prime.md),
[`g.prior()`](http://merliseclyde.github.io/BAS/reference/g.prior.md),
[`hyper.g()`](http://merliseclyde.github.io/BAS/reference/hyper.g.md),
[`hyper.g.n()`](http://merliseclyde.github.io/BAS/reference/hyper.g.n.md),
[`robust()`](http://merliseclyde.github.io/BAS/reference/robust.md),
[`tCCH()`](http://merliseclyde.github.io/BAS/reference/tCCH.md),
[`testBF.prior()`](http://merliseclyde.github.io/BAS/reference/testBF.prior.md)

## Author

Merlise A Clyde

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
