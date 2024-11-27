# Notes to CRAN

## Submission reason 

Bugs: 

- Corrected Additional Issues: C23 (Checks of compiling
C code in C23 mode) on CRAN check page for compiling C code with
clang19 reported 11/16/24. Removed legacy code that was causing the error and if left
would trigger archival of package on CRAN. Github Issue #89

- Fixed Issue #87: Prior inclusions probabilities appear incorrect for a Bernoulli(.2) prior when always including one other predictor bug

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



