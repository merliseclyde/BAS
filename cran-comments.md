# Notes to CRAN

## Submission reason 

This submission addresses CRAN Notes and a bug fix.

CRAN NOTES: 
Version: 1.7.5
Check: compiled code
Result: NOTE 
  File ‘BAS/libs/BAS.so’:
    Found non-API call to R: ‘SETLENGTH’

Replaced use of  non-API call to `SETLENGTH` for over allocated vectors and objects in `C` code for 
sampling without enumeration - Issue #82

Bugs: 

- Fixed Issue #96: where `R2` was incorrect for models where the number of columns in the design matrix exceeded `n`, but
the model was full rank. 

- Fixed Issue #97: where the truncated poisson and truncated power priors did not account for the number of 
models of a given size.

## Test environments

- local OS X install, R 4.5.2 (arm64)
- ubuntu  (github actions CI), R-release R-devel R-oldrelease
- windows (github actions CI), R-release;
- win-builder (r-release, r-devel)

## R CMD check results for this submission

* Mac, Windows, Ubuntu
 0 error | 0 warnings | 0 notes


## Reverse Dependencies

- ginormal
- EMJMCMC
- PEPBVS
- FBMS

## revdepcheck results

## revdepcheck results

We checked 4 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


