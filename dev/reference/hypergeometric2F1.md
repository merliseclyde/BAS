# Gaussian hypergeometric2F1 function

Compute the Gaussian Hypergeometric2F1 function: 2F1(a,b,c,z) =
Gamma(b-c) Int_0^1 t^(b-1) (1 - t)^(c -b -1) (1 - t z)^(-a) dt

## Usage

``` r
hypergeometric2F1(a, b, c, z, method = "Cephes", log = TRUE)
```

## Arguments

- a:

  arbitrary

- b:

  Must be greater 0

- c:

  Must be greater than b if \|z\| \< 1, and c \> b + a if z = 1

- z:

  \|z\| \<= 1

- method:

  The default is to use the Cephes library routine. This sometimes is
  unstable for large a or z near one returning Inf or negative values.
  In this case, try method="Laplace", which use a Laplace approximation
  for tau = exp(t/(1-t)).

- log:

  if TRUE, return log(2F1)

## Value

if log=T returns the log of the 2F1 function; otherwise the 2F1
function.

## Details

The default is to use the routine hyp2f1.c from the Cephes library. If
that return a negative value or Inf, one should try method="Laplace"
which is based on the Laplace approximation as described in Liang et al
JASA 2008. This is used in the hyper-g prior to calculate marginal
likelihoods.

## References

Cephes library hyp2f1.c

Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2005)
Mixtures of g-priors for Bayesian Variable Selection. Journal of the
American Statistical Association. 103:410-423.  
[doi:10.1198/016214507000001337](https://doi.org/10.1198/016214507000001337)

## See also

Other special functions:
[`hypergeometric1F1()`](http://merliseclyde.github.io/BAS/dev/reference/hypergeometric1F1.md),
[`phi1()`](http://merliseclyde.github.io/BAS/dev/reference/phi1.md),
[`trCCH()`](http://merliseclyde.github.io/BAS/dev/reference/trCCH.md)

## Author

Merlise Clyde (<clyde@duke.edu>)

## Examples

``` r
hypergeometric2F1(12, 1, 2, .65)
#> [1] 9.580921
```
