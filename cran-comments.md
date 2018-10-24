## Test environments

* local OS X install, R 3.5.0 R 3.5.1
* ubuntu 14.04.5 (on travis-ci), R 3.5.1 and R-devel
* build-win (devel and release)

## R CMD check results
There were no ERRORs or NOTES.  Warning on win-builder about vignette requiring UTF-8 although locale is ASCII.

## Reverse Dependencies

None


## Comments

added extensive suite of unit tests using testthat in tests directory in order to get code coverage up to 90+ percent.  (part of CII Best Practices)

## Features

* Included an option `pivot=TRUE` in `bas.lm` to fit the models using a pivoted Cholesky decomposition to allow models that are rank-deficient.  [Enhancment #24](https://github.com/merliseclyde/BAS/issues/24) and [Bug #21](https://github.com/merliseclyde/BAS/issues/21).  Currently coefficients that are not-estimable are set to zero so that `predict` and other methods will work as before.  With more testing and timing this may become the default; otherwise the default method without pivoting issues a warning if log marginals are `NA`.  The vector `rank` is added to the output (see documenation for `bas.lm`) and the degrees of freedom methods that assume a uniform prior for obtaining estimates (AIC and BIC) are adjusted to use `rank` rather than `size`.  

* Added option `force.heredity=TRUE`to force lower order terms to be included if higher order terms are present (hierarchical constraint) for `method='MCMC'` and `method='BAS'` with `bas.lm` and `bas.glm`.  Updated Vignette to illustrate. [enhancement #19](https://github.com/merliseclyde/BAS/issues/19).  Checks to see if  _parents_ are included using `include.always` pass  [issue #26](https://github.com/merliseclyde/BAS/issues/26).

* Added option `drop.always.included` to `image.bas` so that variables that are always included may be excluded from the image. By default all are shown [enhancement #23](https://github.com/merliseclyde/BAS/issues/23)

* Added option `drop.always.included` and `subset` to `plot.bas` so that variables that are always included may be excluded from the plot showing the marginal posterior inclusion probabilities (`which=4`). By default all are shown [enhancement #23](https://github.com/merliseclyde/BAS/issues/23)

* update `fitted.bas` to use predict so that code covers both GLM and LM cases with `type='link'` or `type='response'`

* Updates to package for [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/2055/badge)](https://bestpractices.coreinfrastructure.org/projects/2055)  certification

* Added [Code Coverage](https://codecov.io/gh/merliseclyde/BAS) support and more extensive tests using `test_that`.

## Bugs

* fixed [issue #36](https://github.com/merliseclyde/BAS/issues/36) Errors in prior = "ZS-null" when R2 is not finite or out of range due to model being not full rank. Change  in `gexpectations` function in file `bayesreg.c`

* fixed  [issue #35](https://github.com/merliseclyde/BAS/issues/35) for `method="MCMC+BAS"` in `bas.glm`  in `glm_mcmcbas.c` when no values are provided for `MCMC.iterations` or `n.models` and defaults are used.  Added unit test in `test-bas-glm.R`

* fixed  [issue #34](https://github.com/merliseclyde/BAS/issues/34) for `bas.glm` where variables in `include.always` had marginal inclusion probabilities that were incorrect.  Added unit test in `test-bas-glm.R`

* fixed [issue #33](https://github.com/merliseclyde/BAS/issues/33)  for Jeffreys prior where marginal inclusion probabilities were not renomalized after dropping intercept model 

* fixed [issue #32](https://github.com/merliseclyde/BAS/issues/32) 
to allow vectorization for `phi1` function in R/cch.R
and added unit test to "tests/testthat/test-special-functions.R"

* fixed  [issue #31](https://github.com/merliseclyde/BAS/issues/30) to coerce `g` to be a REAL for `g.prior` prior and `IC.prior` in `bas.glm`; added unit-test "tests/testthat/test-bas-glm.R"

* fixed  [issue #30](https://github.com/merliseclyde/BAS/issues/30) added n as hyperparameter if NULL and coerced to be a REAL for `intrinsic` prior in `bas.glm`; added unit-test

* fixed  [issue #29](https://github.com/merliseclyde/BAS/issues/29) added n as hyperparameter if NULL and coerced to be a REAL for `beta.prime` prior in `bas.glm`; added unit-test

* fixed  [issue #28](https://github.com/merliseclyde/BAS/issues/28)  fixed length of MCMC estimates of marginal inclusion probabilities; added unit-test


* fixed  [issue #27](https://github.com/merliseclyde/BAS/issues/27) where expected shrinkage with the JZS prior was greater than 1.  Added unit test.

* fixed output `include.always` to include the intercept [issue #26](https://github.com/merliseclyde/BAS/issues/26) always so that `drop.always.included = TRUE` drops the intercept and any other variables that are forced in.    `include.always` and `force.heredity=TRUE` can now be used together with `method="BAS"`.

* added warning if marginal likelihoods/posterior probabilities are NA with default model fitting method with suggestion that models be rerun with `pivot = TRUE`.  This uses a modified Cholesky decomposition with pivoting so that if the model is rank deficient or nearly singular the dimensionality is reduced.  [Bug #21](https://github.com/merliseclyde/BAS/issues/21).   

* corrected count for first model with `method='MCMC'` which lead to potential model with 0 probabiliy and errors in `image`.

* coerced predicted values to be a vector under BMA (was a matrix)

* fixed `size` with using `method=deterministic` in `bas.glm` (was not updated) 

## Other

*  suppress `warning` when sampling probabilities are 1 or 0 and the number of models is decremented  
[Issue #25](https://github.com/merliseclyde/BAS/issues/25)

* changed `force.heredity.bas` to  renormalize the prior probabilities  rather than to use a new prior probability based on heredity constraints.  For future,  add new priors for models based on heredity.  See comment on  [issue #26](https://github.com/merliseclyde/BAS/issues/26).

* Changed License to GPL 3.0
