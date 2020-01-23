# BAS 1.5.5  Comments to CRAN

resubmission to fix additional errors found during Version 1.5.4  (1-20-2020)

Flavor: r-devel-linux-x86_64-fedora-clang
Check: re-building of vignette outputs 
Result: WARN 
    Error(s) in re-building vignettes:
    --- re-building ‘BAS-vignette.Rmd’ using rmarkdown
    Quitting from lines 463-466 (BAS-vignette.Rmd) 
    Error: processing vignette 'BAS-vignette.Rmd' failed with diagnostics:
    cannot open the connection to 'https://stat.duke.edu/sites/stat.duke.edu/files/climate.dat'
    --- failed re-building ‘BAS-vignette.Rmd’
    

Solaris:
Check: tests 
Result: ERROR 
     Running ‘testthat.R’ [31s/35s]
    Running the tests in ‘tests/testthat.R’ failed.

## Test environments

* local OS X install, R 3.6.2
* debian/gcc-9 fortran-9/valgrind R-devel (4.0.0)
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


## valgrind report (debian/gcc-9/R-devel)
## 
R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes --log-fd=1 --log-file=BAS.Rcheck/BAS-valgrind.txt" --vanilla < BAS.Rcheck/BAS-Ex.R
 
 =1802== Memcheck, a memory error detector
==1802== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==1802== Using Valgrind-3.15.0 and LibVEX; rerun with -h for copyright info
==1802== Command: /usr/local/lib/R/bin/exec/R --vanilla
==1802== Parent PID: 1
==1802==
==1802==
==1802== HEAP SUMMARY:
==1802==     in use at exit: 148,567,369 bytes in 51,660 blocks
==1802==   total heap usage: 473,272 allocs, 421,612 frees, 1,187,563,817 bytes allocated
==1802==
==1802== LEAK SUMMARY:
==1802==    definitely lost: 0 bytes in 0 blocks
==1802==    indirectly lost: 0 bytes in 0 blocks
==1802==      possibly lost: 0 bytes in 0 blocks
==1802==    still reachable: 148,567,369 bytes in 51,660 blocks
==1802==                       of which reachable via heuristic:
==1802==                         newarray           : 4,264 bytes in 1 blocks
==1802==         suppressed: 0 bytes in 0 blocks
==1802== Reachable blocks (those to which a pointer was found) are not shown.
==1802== To see them, rerun with: --leak-check=full --show-leak-kinds=all
==1802==
==1802== For lists of detected and suppressed errors, rerun with: -s
==1802== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

## Reverse Dependencies

 
None

