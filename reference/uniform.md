# Uniform Prior Distribution for Models

Creates an object representing the prior distribution on models for BAS.

## Usage

``` r
uniform()
```

## Value

returns an object of class "prior", with the family name Uniform.

## Details

The Uniform prior distribution is a commonly used prior in BMA, and is a
special case of the independent Bernoulli prior with probs=.5. The
implied prior distribution on model size is binomial(p, .5).

## See also

[`bas.lm`](http://merliseclyde.github.io/BAS/reference/bas.lm.md),
[`beta.binomial`](http://merliseclyde.github.io/BAS/reference/beta.binomial.md),[`Bernoulli`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),

Other priors modelpriors:
[`Bernoulli()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),
[`Bernoulli.heredity()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.heredity.md),
[`beta.binomial()`](http://merliseclyde.github.io/BAS/reference/beta.binomial.md),
[`tr.beta.binomial()`](http://merliseclyde.github.io/BAS/reference/tr.beta.binomial.md),
[`tr.poisson()`](http://merliseclyde.github.io/BAS/reference/tr.poisson.md),
[`tr.power.prior()`](http://merliseclyde.github.io/BAS/reference/tr.power.prior.md)

## Author

Merlise Clyde

## Examples

``` r
uniform()
#> $family
#> [1] "Uniform"
#> 
#> $hyper.parameters
#> [1] 0.5
#> 
#> attr(,"class")
#> [1] "prior"
```
