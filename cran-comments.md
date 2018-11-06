#  BAS 1.5.4 Comments to CRAN

Fixed remaining errors ASAP identified on cran checks   https://cran.r-project.org/web/checks/check_results_BAS.html caught by unit tests on Solaris.  Error in test suite on solaris and fedora/clang depended on state of the random seed and was not immediatly reporducible across platforms, but once identified I was able to identify the issue within the C code and fix. 


## Test environments

* local OS X install, R 3.5.1  (clang / gfortran 6.1)
* local fedora 25  R 3.5.1 (clang)
* ubuntu 14.04.5 (on travis-ci), R 3.5.1 and R-devel  (gcc)
* win-builder (devel and release)
* rhub:  check_with_valgrind (debian+gcc), 
         check_for_cran(platforms = "fedora-clang-devel", "solaris-x86-patched")

## R CMD check results
There were no ERRORs or WARNINGS. 

NOTES:
N  checking CRAN incoming feasibility (47.7s)
   Maintainer: ‘Merlise Clyde <clyde@duke.edu>’
   
   Days since last update: 4  (explained above)

## Reverse Dependencies

None



## News Entry for  BAS 1.5.4
##  Bug Fixes
	* rounding issues with clang on fedora and solaris with
	`force.heredity = TRUE` lead to sampling continuing under
	`method='BAS'` and duplicate models so that normalized posterior
	probabilities were incorrect.
	[issue #38](https://github.com/merliseclyde/BAS/issues/38)

  * FORTRAN errors when data has zero rows 
  (issue #37)[https://github.com/merliseclyde/BAS/issues/37]
  add check and test for n == 0 due to subsetting input data
  


