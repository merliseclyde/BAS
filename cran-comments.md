# BAS 1.6.1 Comments to CRAN

## Test environments

* local OS X install, R 4.1.3
* ubuntu  (on travis-ci, github), R-release R-devel R-oldrelease
* rhub Fedora, Debian, Ubunto, Windows (r-devel)
* win-builder (r-release, r-devel)

## R CMD check results for this submission

* Platform Mac OSX R-release, clang, gfortan 0 error | 0 warnings | 0 notes 

* Platform:   Fedora Linux, R-devel, clang, gfortran  0 error | 0 warnings | 0 notes  
  
* Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN  0 error | 0 warnings | 0 notes  

* Platform:   Windows via rhub (r-devel)  0 errors | 0 warnings  | 1 notes 

── BAS 1.6.1: NOTE

  Build ID:   BAS_1.6.1.tar.gz-26edea4357fb43adb8e7436d72198985
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  22m 21.5s ago
  Build time: 5m 37.5s

> checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✓ | 0 warnings ✓ | 1 note x

NOTE to CRAN:  Note occurs only on Windows r-devel on rhub

* Platform Windows via win-builder (r-devel, r-release)  0 errors | 0 warnings  | 1 notes 

checking CRAN incoming feasibility ... NOTE
Maintainer: 'Merlise Clyde <clyde@duke.edu>'

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1214/088342304000000035
    From: man/BAS.Rd
          man/bas.lm.Rd
    Status: 500
    Message: Internal Server Error
  URL: https://doi.org/10.1214/ss/1009212519
    From: man/bas.lm.Rd
    Status: 500
    Message: Internal Server Error

NOTE to CRAN:  url is valid,  but encounters time-out during check (passes check on all other platforms, and no issues with other urls in package)


## Reverse Dependencies

 
None

