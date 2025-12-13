# eplogprob - Compute approximate marginal inclusion probabilities from pvalues

`eplogprob` calculates approximate marginal posterior inclusion
probabilities from p-values computed from a linear model using a lower
bound approximation to Bayes factors. Used to obtain initial inclusion
probabilities for sampling using Bayesian Adaptive Sampling `bas.lm`

## Usage

``` r
eplogprob(lm.obj, thresh = 0.5, max = 0.99, int = TRUE)
```

## Arguments

- lm.obj:

  a linear model object

- thresh:

  the value of the inclusion probability when if the p-value \>
  1/exp(1), where the lower bound approximation is not valid.

- max:

  maximum value of the inclusion probability; used for the `bas.lm`
  function to keep initial inclusion probabilities away from 1.

- int:

  If the Intercept is included in the linear model, set the marginal
  inclusion probability corresponding to the intercept to 1

## Value

`eplogprob` returns a vector of marginal posterior inclusion
probabilities for each of the variables in the linear model. If int =
TRUE, then the inclusion probability for the intercept is set to 1. If
the model is not full rank, variables that are linearly dependent base
on the QR factorization will have NA for their p-values. In bas.lm,
where the probabilities are used for sampling, the inclusion probability
is set to 0.

## Details

Sellke, Bayarri and Berger (2001) provide a simple calibration of
p-values

BF(p) = -e p log(p)

which provide a lower bound to a Bayes factor for comparing H0: beta = 0
versus H1: beta not equal to 0, when the p-value p is less than 1/e.
Using equal prior odds on the hypotheses H0 and H1, the approximate
marginal posterior inclusion probability

p(beta != 0 \| data ) = 1/(1 + BF(p))

When p \> 1/e, we set the marginal inclusion probability to 0.5 or the
value given by `thresh`.

## References

Sellke, Thomas, Bayarri, M. J., and Berger, James O. (2001),
“Calibration of p-values for testing precise null hypotheses”, The
American Statistician, 55, 62-71.

## See also

[`bas`](http://merliseclyde.github.io/BAS/reference/bas.lm.md)

## Author

Merlise Clyde <clyde@stat.duke.edu>

## Examples

``` r
library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
eplogprob(lm(y ~ ., data=UScrime))
#> (Intercept)           M          So          Ed         Po1         Po2 
#>   1.0000000   0.9480823   0.5000000   0.9883193   0.5045770   0.5000000 
#>          LF         M.F         Pop          NW          U1          U2 
#>   0.5000000   0.5304659   0.5807429   0.8046414   0.5000000   0.6858642 
#>         GDP        Ineq        Prob        Time 
#>   0.6069289   0.9900000   0.9412475   0.5782754 

```
