# Notes to CRAN

## Submission reason 

BAS 1.7.2 which was submitted on 9/16 to 
address the "Strict" and "Notes" issues for R-Devel on the [Checks page at CRAN] (https://cran.r-project.org/web/checks/check_results_BAS.html):

- Removed legacy definitions of ‘PI’ and ‘Free’ and replaced with‘M_PI’ and ‘R_Free’ to comply with ‘STRICT_R_HEADERS’ (issue #81) to prevent package  removal after 9/23/2024

- Removed non-API calls to `SETLENGTH` (issue #82)

The latter change inadvertently is leading to stack imbalances that are not caught in R CMD checks on R-Devel or R-Release with M1Mac but appear in interactive sessions or under checks with valgrind. 

To eliminate the stack imbalance the src code has reverted back to using `SETLENGTH` to truncate over-allocated vectors prior to assigning to the returned SEXP and UNPROTECTING the memory before function exit. (As the number of entries is random, it is not possible to pre-determine the vector size, but just an upper limit.) This code has been in the package since 2018 and further work is still needed to update without creating other memory issues in large problems in practice.

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

This is under development and I ask that CRAN allow this update despite the NOTE to allow the package to be retained on CRAN. The non-API calls identified by `Strict` have been fixed.

## Reverse Dependencies

- ginormal
- EMJMCMC

### revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


