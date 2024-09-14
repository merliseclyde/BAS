# Notes to CRAN

## Submission reason 

- Removed legacy definitions of ‘PI’ and ‘Free’ and replaced with‘M_PI’ and ‘R_Free’ to comply with ‘STRICT_R_HEADERS’ so that package not removed 9/23/2024

## Test environments

- r-devel with valgrind via rhub github actions
- local OS X install, R 4.4.1 (x86, arm64)
- ubuntu  (github actions CI), R-release R-devel R-oldrelease
- win-builder (r-release, r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)
- R-hub linux-x86_64-rocker-gcc-san (r-devel)


## R CMD check results for this submission

* Mmac, Windows, Ubunto, Debian
 0 error | 0 warnings | 0 notes


## Reverse Dependencies

- ginormal
- EMJMCMC

### revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


