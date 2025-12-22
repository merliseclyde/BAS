# Notes to CRAN

## Submission reason 

* This submission addresses CRAN Additional Issues under 

rchk: 
Package BAS version 2.0.0
Package built using 89195/R 4.6.0; x86_64-pc-linux-gnu; 2025-12-19 19:28:42 UTC; unix   
Checked with rchk version 35618ebbccf3cd0b45a3530e6303970a22a9056b LLVM version 14.0.6
More information at https://github.com/kalibera/cran-checks/blob/master/rchk/PROTECT.md
For rchk in docker image see https://github.com/kalibera/rchk/blob/master/doc/DOCKER.md

Function compute_margprobs_Bayes_BAS_MCMC
  [UP] unprotected variable samplemargs while calling allocating function Rprintf BAS/src/model_probabilities.c:134
  [UP] unprotected variable samplemargs while calling allocating function Rprintf BAS/src/model_probabilities.c:149

valgrind:  
unintialized values leading to conditional
jumps or moves depending on unitialized values

valgrind time: 1.5 hours so may time-out on precheck servers

* Update is necessary before 1/15/2026 to maintain on CRAN


## Test environments

- ubuntu R-devel with valgrind via docker 
- local OS X install, R 4.5.2 (arm64)
- win-builder (r-release, r-devel)

No issues identified via the docker valgrind checks and running of examples and unit tests per WRE usage of Valgrind

## R CMD check results for this submission

* Mac, Windows, Ubuntu
 0 error | 0 warnings | 1 notes

Days since last update: 6

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


