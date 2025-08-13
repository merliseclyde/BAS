# Notes to CRAN

## Submission reason 

This submission addresses CRAN Notes and a bug fix.

CRAN Notes: replaced use of  non-API call to `SETLENGTH` for over allocated vectors and objects in `C` code for sampling without enumeration.

Bugs: 

- Fixed Issue #96: where `R2` was incorrect for models where the number of columns in the design matrix exceeded `n`, but
the model was full rank. 

## Test environments

- r-devel with valgrind via rhub github actions
- local OS X install, R 4.4.2 (x86, arm64)
- ubuntu  (github actions CI), R-release R-devel R-oldrelease
- windows (github actions CI), R-release; 
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
I have made some progress in removing the note but it is not yet complete.

The C23 error under Additional Issues has been fixed with removal of legacy code.

## Reverse Dependencies

- ginormal
- EMJMCMC
- PEPBVS

## revdepcheck results

We checked 3 reverse dependencies, comparing R CMD check results across CRAN 
and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages



