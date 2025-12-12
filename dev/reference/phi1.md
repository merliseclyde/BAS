# Compound Confluent hypergeometric function of two variables

Compute the Confluent Hypergeometric function of two variables, also
know as a Horn hypergeometric function or Humbert's hypergeometric used
in Gordy (1998) with integral representation:

## Usage

``` r
phi1(a, b, c, x, y, log = FALSE)
```

## Arguments

- a:

  a \> 0

- b:

  arbitrary

- c:

  c \> 0

- x:

  x \> 0

- y:

  y \> 0

- log:

  logical indicating whether to return phi1 on the log scale

## Details

phi_1(a,b,c,x,y) = \[(Gamma(c)/Gamma(a) Gamma(a-c))\] Int_0^1 t^(a-1)
(1 - t)^(c-a-1) (1 - yt)^(-b) exp(x t) dt
<https://en.wikipedia.org/wiki/Humbert_series> Note that Gordy's
arguments for x and y are reversed in the reference above.

The original \`phi1\` function in \`BAS\` was based on \`C\` code
provided by Gordy. This function returns NA's when x is greater than
\`log(.Machine\$double.xmax)/2\`. A more stable method for calculating
the \`phi1\` function using R's \`integrate\` was suggested by Daniel
Heemann and is now an option whenever \$x\$ is too large. For
calculating Bayes factors that use the \`phi1\` function we recommend
using the \`log=TRUE\` option to compute log Bayes factors.

## References

Gordy 1998

## See also

Other special functions:
[`hypergeometric1F1()`](http://merliseclyde.github.io/BAS/dev/reference/hypergeometric1F1.md),
[`hypergeometric2F1()`](http://merliseclyde.github.io/BAS/dev/reference/hypergeometric2F1.md),
[`trCCH()`](http://merliseclyde.github.io/BAS/dev/reference/trCCH.md)

## Author

Merlise Clyde (<clyde@duke.edu>)

Daniel Heemann (<df.heemann@gmail.com>)

## Examples

``` r
# special cases
# phi1(a, b, c, x=0, y) is the same as 2F1(b, a; c, y)
phi1(1, 2, 1.5, 0, 1 / 100, log=FALSE)
#> [1] 1.013495
hypergeometric2F1(2, 1, 1.5, 1 / 100, log = FALSE)
#> [1] 1.013495

# phi1(a,0,c,x,y) is the same as 1F1(a,c,x)
phi1(1, 0, 1.5, 3, 1 / 100)
#> [1] 10.13001
hypergeometric1F1(1, 1.5, 3, log = FALSE)
#> [1] 10.13001

# use direct integration
phi1(1, 2, 1.5, 1000, 0, log=TRUE)
#> [1] 996.4253
```
