# Confluent hypergeometric1F1 function

Compute the Confluent Hypergeometric function: 1F1(a,b,c,t) =
Gamma(b)/(Gamma(b-a)Gamma(a)) Int_0^1 t^(a-1) (1 - t)^(b-a-1) exp(c t)
dt

## Usage

``` r
hypergeometric1F1(a, b, c, laplace = FALSE, log = TRUE)
```

## Arguments

- a:

  arbitrary

- b:

  Must be greater 0

- c:

  arbitrary

- laplace:

  The default is to use the Cephes library; for large a or s this may
  return an NA, Inf or negative values,, in which case you should use
  the Laplace approximation.

- log:

  if TRUE, return log(1F1)

## References

Cephes library hyp1f1.c

## See also

Other special functions:
[`hypergeometric2F1()`](http://merliseclyde.github.io/BAS/reference/hypergeometric2F1.md),
[`phi1()`](http://merliseclyde.github.io/BAS/reference/phi1.md),
[`trCCH()`](http://merliseclyde.github.io/BAS/reference/trCCH.md)

## Author

Merlise Clyde (<clyde@stat.duke.edu>)

## Examples

``` r
hypergeometric1F1(11.14756, 0.5, 0.00175097)
#> [1] 0.03856253

```
