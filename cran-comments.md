# Notes to CRAN

## Submission reason 

* This submission addresses CRAN Additional Issues under valgrind:  unintialized values leading to conditional
jumps or moves depending on unitialized values

* Update is necessary before 1/15/2026 to maintain on CRAN


## Test environments

- ubuntu R-devel with valgrind via docker 
- local OS X install, R 4.5.2 (arm64)
- win-builder (r-release, r-devel)

No issues identified via the docker valgrind checks and running of examples and unit tests per WRE usage of Valgrind

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


