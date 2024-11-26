# Notes to CRAN

## Submission reason 

Bugs: 

- fixed issue #89: ERROR noted under Additional Issues: C23 (Checks of compiling
C code in C23 mode) on CRAN check page for compiling C code with
clang19 reported 11/16/24. Removed legacy code that was causing the error.

- Fixed Issue #87: Prior inclusions probabilities appear incorrect for a Bernoulli(.2) prior when always including one other predictor bug

## Test environments

- r-devel with valgrind via rhub github actions
- local OS X install, R 4.4.1 (x86, arm64)
- ubuntu  (github actions CI), R-release R-devel R-oldrelease
- win-builder (r-release, r-devel)

## R CMD check results for this submission

* Mac, Windows, Ubunto
 0 error | 0 warnings | 1 notes

* checking compiled code ... NOTE
File 'BAS/libs/x64/BAS.dll':
  Found non-API call to R: 'SETLENGTH'

Compiled code should not call non-API entry points in R.

This is under development and I ask that CRAN allow this update 
despite the NOTE to allow the package to be retained on CRAN. 
The non-API calls identified by `Strict` have been fixed.

## Reverse Dependencies

- ginormal
- EMJMCMC

### revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


