# Notes to CRAN

## Submission reason

* Subbmission addresses NOTE on CRAN Check Restults
Version: 2.0.2
Check: R code for possible problems
Result: NOTE 
  Found calls to structure() using deprecated special names:
    BAS/tests/testthat/test-bas-glm.R (.Label: 3)
  '.Label' should be changed to 'levels'.
Flavors: r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-debian-gcc, r-devel-linux-x86_64-fedora-clang, 
r-devel-linux-x86_64-fedora-gcc, r-devel-windows-x86_64

* additional issues bugs in code for `bas.lm` and `bas.glm` using weights 
were addressed.


## Test environments

- local OS X install, R 4.6.2 (arm64)
- win-builder (r-release, r-devel)
- gcc (Ubuntu 15.2.0-16ubuntu1)  R-devel

No issues identified via checks and running of examples and unit tests

## R CMD check results for this submission

* Mac, Windows, Ubuntu
 0 error | 0 warnings | 1 notes
 
 Maintainer: 'Merlise Clyde <clyde@duke.edu>'


## Reverse Dependencies

- ginormal
- EMJMCMC
- PEPBVS
- FBMS


## revdepcheck results

We checked 4 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


