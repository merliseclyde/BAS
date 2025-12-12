# Fitting Generalized Linear Models and Bayesian marginal likelihood evaluation

A version of glm.fit rewritten in C; also returns marginal likelihoods
for Bayesian model comparison

## Usage

``` r
bayesglm.fit(
  x,
  y,
  weights = rep(1, nobs),
  start = NULL,
  etastart = NULL,
  mustart = NULL,
  offset = rep(0, nobs),
  family = binomial(),
  coefprior = bic.prior(nobs),
  control = glm.control(),
  intercept = TRUE
)
```

## Arguments

- x:

  design matrix

- y:

  response

- weights:

  optional vector of weights to be used in the fitting process. Should
  be NULL or a numeric vector.

- start:

  starting value for coefficients in the linear predictor

- etastart:

  starting values for the linear predictor

- mustart:

  starting values for the vectors of means

- offset:

  a priori known component to be included in the linear predictor

- family:

  a description of the error distribution and link function for
  exponential family; currently only binomial(), poisson(), and Gamma()
  with canonical links are implemented.

- coefprior:

  function specifying prior distribution on coefficients with optional
  hyperparameters leading to marginal likelihood calculations; options
  include
  [`bic.prior()`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md),` `[`aic.prior()`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md),
  and `ic.prior()`

- control:

  a list of parameters that control convergence in the fitting process.
  See the documentation for
  [`glm.control()`](https://rdrr.io/r/stats/glm.control.html)

- intercept:

  should an intercept be included in the null model?

## Value

- coefficients:

  MLEs

- se:

  Standard errors of coefficients based on the sqrt of the diagonal of
  the inverse information matrix

- mu:

  fitted mean

- rank:

  numeric rank of the fitted linear model

- deviance:

  minus twice the log likelihood evaluated at the MLEs

- g:

  value of g in g-priors

- shrinkage:

  shrinkage factor for coefficients in linear predictor

- RegSS:

  quadratic form beta'I(beta)beta used in shrinkage

- logmarglik:

  the log marginal or integrated log likelihood (up to a constant)

## Details

C version of glm-fit. For different prior choices returns, marginal
likelihood of model using a Laplace approximation.

## References

[`glm`](https://rdrr.io/r/stats/glm.html)

## See also

[`bic.prior`](http://merliseclyde.github.io/BAS/dev/reference/IC.prior.md)

## Author

Merlise Clyde translated the
[`glm.fit`](https://rdrr.io/r/stats/glm.html) from R base into C using
the .Call interface

## Examples

``` r
data(Pima.tr, package="MASS")
Y <- as.numeric(Pima.tr$type) - 1
X <- cbind(1, as.matrix(Pima.tr[,1:7]))
out <- bayesglm.fit(X, Y, family=binomial(),coefprior=bic.prior(n=length(Y)))
out$coef
#> [1] -9.773061533  0.103183427  0.032116823 -0.004767542 -0.001916632
#> [6]  0.083623912  1.820410367  0.041183529
out$se
#> [1] 1.770386016 0.064694153 0.006787299 0.018540741 0.022499541 0.042826888
#> [7] 0.665513776 0.022090978
# using built in function
glm(type ~ ., family=binomial(), data=Pima.tr)
#> 
#> Call:  glm(formula = type ~ ., family = binomial(), data = Pima.tr)
#> 
#> Coefficients:
#> (Intercept)        npreg          glu           bp         skin          bmi  
#>   -9.773062     0.103183     0.032117    -0.004768    -0.001917     0.083624  
#>         ped          age  
#>    1.820410     0.041184  
#> 
#> Degrees of Freedom: 199 Total (i.e. Null);  192 Residual
#> Null Deviance:       256.4 
#> Residual Deviance: 178.4     AIC: 194.4
```
