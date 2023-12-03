# BAS 1.6.7 Comments to CRAN

# Notes to CRAN

## Submission reason 


Status  on CRAN check page indicated errors due to `valgrind` -  Conditional jump or move depends on uninitialised value(s)

- Initialized `se`  via `memset` and `disp` in  `fit_glm.c` (issue #72)

- Fixed issue #67 reported on Github


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

* Mmac, Windows,   Ubunto, Debian
 0 error | 0 warnings | 0 notes


## Reverse Dependencies

* ginormal

### revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


