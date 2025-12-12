# Truncated Compound Confluent Hypergeometric function

Compute the Truncated Confluent Hypergeometric function from Li and
Clyde (2018) which is the normalizing constant in the tcch density of
Gordy (1998) with integral representation:

## Usage

``` r
trCCH(a, b, r, s, v, k, log = FALSE)
```

## Arguments

- a:

  a \> 0

- b:

  b \> 0

- r:

  r \>= 0

- s:

  arbitrary

- v:

  0 \< v

- k:

  arbitrary

- log:

  logical indicating whether to return values on the log scale; useful
  for Bayes Factor calculations

## Details

tr.cch(a,b,r,s,v,k) = Int_0^1/v u^(a-1) (1 - vu)^(b -1) (k + (1 -
k)vu)^(-r) exp(-s u) du

This uses a more stable method for calculating the normalizing constant
using R's \`integrate\` function rather than the version in Gordy 1998.
For calculating Bayes factors that use the \`trCCH\` function we
recommend using the \`log=TRUE\` option to compute log Bayes factors.

## References

Gordy 1998 Li & Clyde 2018

## See also

Other special functions:
[`hypergeometric1F1()`](http://merliseclyde.github.io/BAS/dev/reference/hypergeometric1F1.md),
[`hypergeometric2F1()`](http://merliseclyde.github.io/BAS/dev/reference/hypergeometric2F1.md),
[`phi1()`](http://merliseclyde.github.io/BAS/dev/reference/phi1.md)

## Author

Merlise Clyde (<clyde@duke.edu>)

## Examples

``` r
# special cases
# trCCH(a, b, r, s=0, v = 1, k) is the same as
# 2F1(a, r, a + b, 1 - 1/k)*beta(a, b)/k^r

k = 10; a = 1.5; b = 2; r = 2;  
trCCH(a, b, r, s=0, v = 1, k=k) *k^r/beta(a,b)
#> [1] 4.74679
hypergeometric2F1(a, r, a + b, 1 - 1/k, log = FALSE)
#> [1] 4.746772

# trCCH(a,b,0,s,1,1) is the same as 
# beta(a, b) 1F1(a, a + b, -s, log=FALSE)
s = 3; r = 0; v = 1; k = 1
beta(a, b)*hypergeometric1F1(a, a+b, -s, log = FALSE)
#> [1] 0.0923551
trCCH(a, b, r, s, v, k)
#> [1] 0.09235518

# Equivalence with the Phi1 function 
a = 1.5; b = 3; k = 1.25; s = 400;  r = 2;  v = 1; 

phi1(a, r,  a + b, -s, 1 - 1/k,  log=FALSE)*(k^-r)*gamma(a)*gamma(b)/gamma(a+b)
#> [1] 7.04733e-05
trCCH(a,b,r,s,v,k)
#> [1] 7.052186e-05
```
