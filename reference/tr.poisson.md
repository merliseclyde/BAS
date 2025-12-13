# Truncated Poisson Prior Distribution for Models

Creates an object representing the prior distribution on models for BAS
using a truncated Poisson Distribution on the Model Size

## Usage

``` r
tr.poisson(lambda, trunc)
```

## Arguments

- lambda:

  parameter in the Poisson distribution representing expected model size
  with infinite predictors

- trunc:

  parameter that determines truncation in the distribution i.e. P(M;
  lambda, trunc) = 0 if M \> trunc

## Value

returns an object of class "prior", with the family and hyperparameters.

## Details

The Poisson prior distribution on model size is obtained by assigning
each variable inclusion indicator independent Bernoulli distributions
with probability w, and then taking a limit as p goes to infinity and w
goes to zero, such that p\*w converges to lambda. The Truncated version
assigns zero probability to all models of size M \> trunc.

## See also

[`bas.lm`](http://merliseclyde.github.io/BAS/reference/bas.lm.md),
[`Bernoulli`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),[`uniform`](http://merliseclyde.github.io/BAS/reference/uniform.md)

Other priors modelpriors:
[`Bernoulli()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.md),
[`Bernoulli.heredity()`](http://merliseclyde.github.io/BAS/reference/Bernoulli.heredity.md),
[`beta.binomial()`](http://merliseclyde.github.io/BAS/reference/beta.binomial.md),
[`tr.beta.binomial()`](http://merliseclyde.github.io/BAS/reference/tr.beta.binomial.md),
[`tr.power.prior()`](http://merliseclyde.github.io/BAS/reference/tr.power.prior.md),
[`uniform()`](http://merliseclyde.github.io/BAS/reference/uniform.md)

## Author

Merlise Clyde

## Examples

``` r
tr.poisson(10, 50)
#> $family
#> [1] "Trunc-Poisson"
#> 
#> $hyper.parameters
#> [1] 10 50
#> 
#> attr(,"class")
#> [1] "prior"
```
