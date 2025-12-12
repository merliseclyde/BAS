# Coerce a BAS list object into a matrix.

Models, coefficients, and standard errors in objects of class 'bas' are
represented as a list of lists to reduce storage by omitting the zero
entries. These functions coerce the list object to a matrix and fill in
the zeros to facilitate other computations.

## Usage

``` r
list2matrix.which(x, which.models = NULL)
```

## Arguments

- x:

  a 'bas' object

- which.models:

  a vector of indices use to extract a subset

## Value

a matrix representation of `x$what`, with number of rows equal to the
length of which.models or total number of models and number of columns
`x$n.vars`

## Details

`list2matrix.bas(x, which)` is equivalent to `list2matrix.which(x)`,
however, the latter uses sapply rather than a loop. `list2matrix.which`
and `which.matrix` both coerce `x$which` into a matrix.

## See also

[`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md)

Other as.matrix methods:
[`list2matrix.bas()`](http://merliseclyde.github.io/BAS/dev/reference/list2matrix.md),
[`which.matrix()`](http://merliseclyde.github.io/BAS/dev/reference/which.matrix.md)

## Author

Merlise Clyde <clyde@duke.edu>

## Examples

``` r
data(Hald)
Hald.bic <-  bas.lm(Y ~ ., data=Hald, prior="BIC", initprobs="eplogp")
coef <- list2matrix.bas(Hald.bic, "mle")  # extract all ols coefficients
se <- list2matrix.bas(Hald.bic, "mle.se")
models <- list2matrix.which(Hald.bic)     #matrix of model indicators
models <- which.matrix(Hald.bic$which, Hald.bic$n.vars)     #matrix of model indicators
```
