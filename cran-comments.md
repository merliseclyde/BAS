#  BAS 1.5.3 Comments to CRAN

I was requested to fix errors ASAP identified on cran checks   https://cran.r-project.org/web/checks/check_results_BAS.html with the submission 4 days ago.
Version 1.5.2 included an extensive suite of unit tests. While extensive testing on several platforms prior to submission to CRAN suggested that the package passed CRAN checks,  the unit-tests did catch additional errors on debian/fedora using clang post-submission.  Subsequent testing on fedora with clang (locally) and debian/gcc/valgrind via rhub indicate that the package passes checks (as well as on the previous platforms tested).   Remaining results from valgrind are not linked to lines in the package C code, so I have been unable to track those (similar messages appear on startup with R), so believe that they are false positives.  


## Test environments

* local OS X install, R 3.5.0 R 3.5.1  (clang / gfortran 6.1)
* local fedora R 3.5.1 (clang + valgrind)
* ubuntu 14.04.5 (on travis-ci), R 3.5.1 and R-devel  (gcc)
* win-builder (devel and release)
* rhub:  check_with_valgrind (debian+gcc), 
         check_for_cran(platforms = "fedora-clang-devel")

## R CMD check results
There were no ERRORs or WARNINGS. 

NOTES:
N  checking CRAN incoming feasibility (47.7s)
   Maintainer: ‘Merlise Clyde <clyde@duke.edu>’
   
   Days since last update: 4  (explained above)

## Reverse Dependencies

None



## News Entry for  BAS 1.5.3

Fixed errors identified on cran checks https://cran.r-project.org/web/checks/check_results_BAS.html

* initialize R2_m = 0.0 in lm_mcmcbas.c  (lead to NA's with clang on debian and fedora )

* switch to default of `pivot = TRUE` in `bas.lm`, adding `tol` as an argument to control tolerance in `cholregpovot` for improved stability across platforms with singular or nearly singular designs.

* valgrind messages: Conditional jump or move depends on uninitialised value(s). Initialize vectors allocated via R_alloc in lm_deterministic.c and glm_deterministic.c.
