#  BAS 1.5.4 Comments to CRAN
#  
I was requested to fix valgrind errors identified on cran checks   https://cran.r-project.org/web/checks/check_results_BAS.html. 
Subsequent testing on fedora with clang (locally) and debian/gcc/valgrind via rhub indicate that the package passes checks.  

## Test environments

* local OS X install, R 3.6.2
* local fedora R 3.6.1 (clang + valgrind)
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)
* R-hub windows-x86_64-devel (r-devel)
* R-hub ubuntu-gcc-release (r-release)
* R-hub fedora-clang-devel (r-devel)
* R-hub linux-x86_64-rocker-gcc-san (r-devel)
* R-hub solaris-x86-patched (r-patched)
* R-hub debian-gcc-release (r-release) (with valgrind)

## R CMD check results

No Notes, Warnings or Errors on platforms except 

* On fedora-clang-devel (r-devel):  0 error | 0 warnings | 2 notes  (False positives)

 checking CRAN incoming feasibility ...NB: need Internet access to use CRAN incoming checks
  NOTE
  Maintainer: ‘Merlise Clyde <clyde@duke.edu>’
  
  Possibly mis-spelled words in DESCRIPTION:
    BAS (32:30)
    DMS (42:28)
    Ghosh (40:12)
    Liang (25:30)
    Littman (40:22)
    Siow (24:34)
    Zellner (24:26)
    Zellner's (23:10)
    al (25:39)
    et (25:36)

* On Solaris: 1 error | 0 warnings | 0 notes 

Error: testthat unit tests failed Execution halted

This occurred with the previous release and significant debugging has not been able to track down the issue although likely machine precision difference with 32 versus 64 bit precision.  Due to wide spread use of the package we request that it remain on CRAN.

## Reverse Dependencies

 
  None
