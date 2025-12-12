# Images of models used in Bayesian model averaging

Creates an image of the models selected using
[`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md).

## Usage

``` r
# S3 method for class 'bas'
image(
  x,
  top.models = 20,
  intensity = TRUE,
  prob = TRUE,
  log = TRUE,
  rotate = TRUE,
  color = "rainbow",
  subset = NULL,
  drop.always.included = FALSE,
  offset = 0.75,
  digits = 3,
  vlas = 2,
  plas = 0,
  rlas = 0,
  ...
)
```

## Arguments

- x:

  A BMA object of type 'bas' created by BAS

- top.models:

  Number of the top ranked models to plot

- intensity:

  Logical variable, when TRUE image intensity is proportional to the
  probability or log(probability) of the model, when FALSE, intensity is
  binary indicating just presence (light) or absence (dark) of a
  variable.

- prob:

  Logical variable for whether the area in the image for each model
  should be proportional to the posterior probability (or log
  probability) of the model (TRUE) or with equal area (FALSE).

- log:

  Logical variable indicating whether the intensities should be based on
  log posterior odds (TRUE) or posterior probabilities (FALSE). The log
  of the posterior odds is for comparing the each model to the worst
  model in the top.models.

- rotate:

  Should the image of models be rotated so that models are on the y-axis
  and variables are on the x-axis (TRUE)

- color:

  The color scheme for image intensities. The value "rainbow" uses the
  rainbow palette. The value "blackandwhite" produces a black and white
  image (greyscale image)

- subset:

  indices of variables to include/exclude in plot

- drop.always.included:

  logical variable to drop variables that are always forced into the
  model. FALSE by default.

- offset:

  numeric value to add to intensity

- digits:

  number of digits in posterior probabilities to keep

- vlas:

  las parameter for placing variable names; see par

- plas:

  las parameter for posterior probability axis

- rlas:

  las parameter for model ranks

- ...:

  Other parameters to be passed to the `image` and `axis` functions.

## Details

Creates an image of the model space sampled using
[`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md). If a
subset of the top models are plotted, then probabilities are
renormalized over the subset.

## Note

Suggestion to allow area of models be proportional to posterior
probability due to Thomas Lumley

## References

Clyde, M. (1999) Bayesian Model Averaging and Model Search Strategies
(with discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid,
J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages
157-185.

## See also

[`bas`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md)

Other bas methods:
[`BAS`](http://merliseclyde.github.io/BAS/dev/reference/BAS.md),
[`bas.lm()`](http://merliseclyde.github.io/BAS/dev/reference/bas.lm.md),
[`coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/coef.md),
[`confint.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.coef.md),
[`confint.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/confint.pred.md),
[`diagnostics()`](http://merliseclyde.github.io/BAS/dev/reference/diagnostics.md),
[`fitted.bas()`](http://merliseclyde.github.io/BAS/dev/reference/fitted.md),
[`force.heredity.bas()`](http://merliseclyde.github.io/BAS/dev/reference/force.heredity.bas.md),
[`plot.confint.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.confint.md),
[`predict.bas()`](http://merliseclyde.github.io/BAS/dev/reference/predict.bas.md),
[`predict.basglm()`](http://merliseclyde.github.io/BAS/dev/reference/predict.basglm.md),
[`summary.bas()`](http://merliseclyde.github.io/BAS/dev/reference/summary.md),
[`update.bas()`](http://merliseclyde.github.io/BAS/dev/reference/update.md),
[`variable.names.pred.bas()`](http://merliseclyde.github.io/BAS/dev/reference/variable.names.pred.bas.md)

Other bas plots:
[`plot.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.md),
[`plot.coef.bas()`](http://merliseclyde.github.io/BAS/dev/reference/plot.coef.md)

## Author

Merlise Clyde <clyde@stat.duke.edu>

## Examples

``` r
require(graphics)
data("Hald")
hald.ZSprior <- bas.lm(Y ~ ., data = Hald, prior = "ZS-null")
image(hald.ZSprior, drop.always.included = TRUE) # drop the intercept
```
