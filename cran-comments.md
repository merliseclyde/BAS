# BAS 1.6.4 Comments to CRAN

# Submission to address Error in CRAN checks for M1mac

output does not produce NA's on M1Mac so warning should not be triggered and test updated

## Test environments

- mac-builder r-release  macosx-arm64 
- local OS X install, R 4.2.1 (x86)
- ubuntu  (github actions CI), R-release R-devel R-oldrelease
- win-builder (r-release, r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)
- R-hub linux-x86_64-rocker-gcc-san (r-devel)


## R CMD check results for this submission

* Platform Mac OSX R-release, M1mac (arm64) clang, gfortran  0 error | 0 warnings | 0 notes

* Platform Mac OSX R-release, (x86) clang, gfortan 0 error | 0 warnings | 0 notes 

* Platform:   Windows via  winbuilder (r-release, r-devel)  0 errors | 0 warnings  | 0 notes   

* Platform:   Ubuntu Linux, R-release, GCC  0 error | 0 warnings | 0 notes  

* Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN  0 error | 0 warnings | 0 notes  

* Platform:   Fedora Linux, R-devel, clang, gfortran  0 error | 0 warnings | 1 notes  
 
  On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
 
(spurious Note)


## Reverse Dependencies

 
None

