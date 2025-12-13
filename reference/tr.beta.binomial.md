# Truncated Beta-Binomial Prior Distribution for Models

Creates an object representing the prior distribution on models for BAS
using a truncated Beta-Binomial Distribution on the Model Size

## Usage

``` r
tr.beta.binomial(alpha = 1, beta = 1, trunc)
```

## Arguments

- alpha:

  parameter in the beta prior distribution

- beta:

  parameter in the beta prior distribution

- trunc:

  parameter that determines truncation in the distribution i.e. P(M;
  alpha, beta, trunc) = 0 if M \> trunc.

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
[`tr.poisson()`](http://merliseclyde.github.io/BAS/reference/tr.poisson.md),
[`tr.power.prior()`](http://merliseclyde.github.io/BAS/reference/tr.power.prior.md),
[`uniform()`](http://merliseclyde.github.io/BAS/reference/uniform.md)

## Author

Merlise Clyde

## Examples

``` r
tr.beta.binomial(1, 10, 5)
#> $family
#> [1] "Trunc-Beta-Binomial"
#> 
#> $hyper.parameters
#> [1]  1 10  5
#> 
#> attr(,"class")
#> [1] "prior"
library(MASS)
data(UScrime)
UScrime[, -2] <- log(UScrime[, -2])
crime.bic <- bas.lm(y ~ .,
  data = UScrime, n.models = 2^15, prior = "BIC",
  modelprior = tr.beta.binomial(1, 1, 8),
  initprobs = "eplogp"
)
```
