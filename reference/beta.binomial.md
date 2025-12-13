# Beta-Binomial Prior Distribution for Models

Creates an object representing the prior distribution on models for BAS.

## Usage

``` r
beta.binomial(alpha = 1, beta = 1)
```

## Arguments

- alpha:

  parameter in the beta prior distribution

- beta:

  parameter in the beta prior distribution

## Value

returns an object of class "prior", with the family and hyperparameters.

## Details

The beta-binomial distribution on model size is obtained by assigning
each variable inclusion indicator independent Bernoulli distributions
with probability w, and then giving w a beta(alpha,beta) distribution.
Marginalizing over w leads to the distribution on model size having the
beta-binomial distribution. The default hyperparameters lead to a
uniform distribution over model size.

## See also

[`bas.lm`](http://merliseclyde.github.io/BAS/reference/bas.lm.md),
[`Bernoulli`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),[`uniform`](http://merliseclyde.github.io/BAS/reference/uniform.md)

Other priors modelpriors:
[`Bernoulli()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),
[`Bernoulli.heredity()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.heredity.md),
[`tr.beta.binomial()`](http://merliseclyde.github.io/BAS/reference/tr.beta.binomial.md),
[`tr.poisson()`](http://merliseclyde.github.io/BAS/reference/tr.poisson.md),
[`tr.power.prior()`](http://merliseclyde.github.io/BAS/reference/tr.power.prior.md),
[`uniform()`](http://merliseclyde.github.io/BAS/reference/uniform.md)

## Author

Merlise Clyde

## Examples

``` r
beta.binomial(1, 10) #' @family priors modelpriors
#> $family
#> [1] "Beta-Binomial"
#> 
#> $hyper.parameters
#> [1]  1 10
#> 
#> attr(,"class")
#> [1] "prior"
```
