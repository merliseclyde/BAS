# BAS 1.6.5 Comments to CRAN

# Notes to CRAN

## Submission reason 

Status  on CRAN check page indicated 1 Note and up to 15 Warnings

- fixed NOTE from checkRd `bas.lm.Rd:109`: Lost braces; missing escapes or markup?

- fixed all WARNINGs related to error format mismatch of type in C code


Fixed issues #39, #56, #61, #68, #69 reported on Github
and unit tests added


## Test environments

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


