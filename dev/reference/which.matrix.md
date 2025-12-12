# Coerce a BAS list object of models into a matrix.

This function coerces the list object of models to a matrix and fill in
the zeros to facilitate other computations.

## Usage

``` r
which.matrix(which, n.vars)
```

## Arguments

- which:

  a 'bas' model object `x$which`

- n.vars:

  the total number of predictors, `x$n.vars`

## Value

a matrix representation of `x$which`, with number of rows equal to the
length of which.models or total number of models and number of columns
`x$n.vars`

## Details

`which.matrix` coerces `x$which` into a matrix.

## See also

[`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md)

Other as.matrix methods:
[`list2matrix.bas()`](http://merliseclyde.github.io/BAS/dev/reference/list2matrix.md),
[`list2matrix.which()`](http://merliseclyde.github.io/BAS/dev/reference/list2matrix.which.md)

## Author

Merlise Clyde <clyde@duke.edu>

## Examples

``` r
data(Hald)
Hald.bic <-  bas.lm(Y ~ ., data=Hald, prior="BIC", initprobs="eplogp")
# matrix of model indicators
models <- which.matrix(Hald.bic$which, Hald.bic$n.vars)
```
