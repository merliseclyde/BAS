# Independent Bernoulli Prior Distribution for Models

Creates an object representing the prior distribution on models for BAS.

## Usage

``` r
Bernoulli(probs = 0.5)
```

## Arguments

- probs:

  a scalar or vector of prior inclusion probabilities. If a scalar, the
  values is replicated for all variables ans a 1 is added for the
  intercept. BAS checks to see if the length is equal to the dimension
  of the parameter vector for the full model and adds a 1 to include the
  intercept.

## Value

returns an object of class "prior", with the family and hyperparameters.

## Details

The independent Bernoulli prior distribution is a commonly used prior in
BMA, with the Uniform distribution a special case with probs=.5. If all
indicator variables have a independent Bernoulli distributions with
common probability probs, the distribution on model size binomial(p,
probs) distribution.

## See also

[`bas.lm`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`beta.binomial`](http://merliseclyde.github.io/BAS/dev/reference/beta.binomial.md),[`uniform`](http://merliseclyde.github.io/BAS/dev/reference/uniform.md)` `

Other priors modelpriors:
[`Bernoulli.heredity()`](http://merliseclyde.github.io/BAS/dev/reference/Bernoulli.heredity.md),
[`beta.binomial()`](http://merliseclyde.github.io/BAS/dev/reference/beta.binomial.md),
[`tr.beta.binomial()`](http://merliseclyde.github.io/BAS/dev/reference/tr.beta.binomial.md),
[`tr.poisson()`](http://merliseclyde.github.io/BAS/dev/reference/tr.poisson.md),
[`tr.power.prior()`](http://merliseclyde.github.io/BAS/dev/reference/tr.power.prior.md),
[`uniform()`](http://merliseclyde.github.io/BAS/dev/reference/uniform.md)

## Author

Merlise Clyde

## Examples

``` r
Bernoulli(.9)
#> $family
#> [1] "Bernoulli"
#> 
#> $hyper.parameters
#> [1] 0.9
#> 
#> attr(,"class")
#> [1] "prior"
```
