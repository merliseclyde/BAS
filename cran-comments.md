# BAS 1.6.0  Comments to CRAN

## submission to include  USE_FC_LEN_T for calls to Fortran subroutines with character strings in light of requirement under R 4.2.0


## Test environments

* local OS X install, R 4.1.2
* ubuntu  (on travis-ci, github), R-release R-devel R-oldrelease
* rhub Fedora, Debian, Ubunto, Windows (r-devel)
* win-builder (r-release, r-devel)


## R CMD check results for this submission

* Platform:   Fedora Linux, R-devel, clang, gfortran  0 error | 0 warnings | 0 notes  
  
* Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN  0 error | 0 warnings | 0 notes  

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

