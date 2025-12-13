# eplogprob.marg - Compute approximate marginal inclusion probabilities from pvalues

`eplogprob.marg` calculates approximate marginal posterior inclusion
probabilities from p-values computed from a series of simple linear
regression models using a lower bound approximation to Bayes factors.
Used to order variables and if appropriate obtain initial inclusion
probabilities for sampling using Bayesian Adaptive Sampling `bas.lm`

## Usage

``` r
eplogprob.marg(Y, X, thresh = 0.5, max = 0.99, int = TRUE)
```

## Arguments

- Y:

  response variable

- X:

  design matrix with a column of ones for the intercept

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

`eplogprob.prob` returns a vector of marginal posterior inclusion
probabilities for each of the variables in the linear model. If int =
TRUE, then the inclusion probability for the intercept is set to 1.

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
value given by `thresh`. For the eplogprob.marg the marginal p-values
are obtained using statistics from the p simple linear regressions

P(F \> (n-2) R2/(1 - R2)) where F ~ F(1, n-2) where R2 is the square of
the correlation coefficient between y and X_j.

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
