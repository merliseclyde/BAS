# Truncated Power Prior Distribution for Models

Creates an object representing the prior distribution on models for BAS
using a truncated Distribution on the Model Size where the probability
of gamma = p^-kappa \|gamma\| where gamma is the vector of model
indicators

## Usage

``` r
tr.power.prior(kappa = 2, trunc)
```

## Arguments

- kappa:

  parameter in the prior distribution that controls sparsity

- trunc:

  parameter that determines truncation in the distribution i.e. P(gamma;
  alpha, beta, trunc) = 0 if \|gamma\| \> trunc.

## Value

returns an object of class "prior", with the family and hyperparameters.

## Details

The beta-binomial distribution on model size is obtained by assigning
each variable inclusion indicator independent Bernoulli distributions
with probability w, and then giving w a beta(alpha,beta) distribution.
Marginalizing over w leads to the number of included predictors having a
beta-binomial distribution. The default hyperparameters lead to a
uniform distribution over model size. The Truncated version assigns zero
probability to all models of size \> trunc.

## See also

[`bas.lm`](http://merliseclyde.github.io/BAS/reference/bas.lm.md),
[`Bernoulli`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),[`uniform`](http://merliseclyde.github.io/BAS/reference/uniform.md)

Other priors modelpriors:
[`Bernoulli()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),
[`Bernoulli.heredity()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.heredity.md),
[`beta.binomial()`](http://merliseclyde.github.io/BAS/reference/beta.binomial.md),
[`tr.beta.binomial()`](http://merliseclyde.github.io/BAS/reference/tr.beta.binomial.md),
[`tr.poisson()`](http://merliseclyde.github.io/BAS/reference/tr.poisson.md),
[`uniform()`](http://merliseclyde.github.io/BAS/reference/uniform.md)

## Author

Merlise Clyde

## Examples

``` r
tr.power.prior(2, 8)
#> $family
#> [1] "Trunc-Power-Prior"
#> 
#> $hyper.parameters
#> [1] 2 8
#> 
#> attr(,"class")
#> [1] "prior"
library(MASS)
data(UScrime)
UScrime[, -2] <- log(UScrime[, -2])
crime.bic <- bas.lm(y ~ .,
  data = UScrime, n.models = 2^15, prior = "BIC",
  modelprior = tr.power.prior(2, 8),
  initprobs = "eplogp"
)
```
