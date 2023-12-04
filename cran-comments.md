# BAS 1.7.0 Comments to CRAN

# Notes to CRAN

## Submission reason 

Submission prior to Dec 14 required to maintain package on CRAN.
Status  on CRAN check page  under `Additional Issues` for `valgrind`  showed warnings  `Conditional jump or move depends on uninitialised value(s)`.  

- Initialized vector `se`  via `memset` and `disp = 1.0` in  `fit_glm.c` (issue #72)

- Fixed issue #67 reported on Github and added unit test


## Test environments

- r-devel with valgrind via docker
- mac-builder r-release  macosx-arm64 
- local OS X install, R 4.3.2 (x86)
- ubuntu  (github actions CI), R-release R-devel R-oldrelease
- win-builder (r-release, r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)
- R-hub linux-x86_64-rocker-gcc-san (r-devel)


## R CMD check results for this submission

* Mmac, Windows, Ubunto, Debian
 0 error | 0 warnings | 1 notes

Note: Days since last update: 5

(resubmission requested by CRAN)

## Reverse Dependencies

* ginormal

### revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


