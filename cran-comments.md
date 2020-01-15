#  BAS 1.5.4 Comments to CRAN  (resubmission)
 
This submission addresses additional issues found on the *last* version released on CRAN:
 LTO <https://urldefense.proofpoint.com/v2/url?u=https-3A__www.stats.ox.ac.uk_pub_bdr_LTO_BAS.out&d=DwICAw&c=imBPVzF25OnBgGmVOlcsiEgHoG1i6YHLR0Sj_gZ4adc&r=NOkxkvdFOOffXzeTY2kgZQ&m=52W89JkPbOpG86B0ILRzgnU71_hsgNFiQNKrlb5-Nt0&s=r-lKhBavxLKaoTr_9AlOwdBxfkFOR0sTSi-oR-WKzxU&e= >
 valgrind <https://urldefense.proofpoint.com/v2/url?u=https-3A__www.stats.ox.ac.uk_pub_bdr_memtests_valgrind_BAS&d=DwICAw&c=imBPVzF25OnBgGmVOlcsiEgHoG1i6YHLR0Sj_gZ4adc&r=NOkxkvdFOOffXzeTY2kgZQ&m=52W89JkPbOpG86B0ILRzgnU71_hsgNFiQNKrlb5-Nt0&s=xZMf4TAeqaSYWV8VES3z6OGZLkp320XviDiKS8Wj8Yg&e= >

Testing with debian/gcc-9/gfortran-9/Valgrind-3.15.0  via rocker/r-devel indicate that the package passes above checks for valgrind (see below) and LTO.


The ERRORS in the *last* submission for Solaris remain (errors were allowed for last submission).

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

This occurred with the previous release and significant debugging has not been able to track down the issue although likely machine precision difference with 32 versus 64 bit precision.  Due to wide spread use of the package we request that it remain on CRAN.

## Other Issues: valgrind report (debian/gcc-9/R-devel)
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

## NEWS for BAS 1.5.4

## Features

* Modified prior probabilities to adjust for the number of variables always
included when using include.always.  [Pull request #41](https://github.com/merliseclyde/BAS/pull/41) by Don van de Bergh.  [Issue #40](https://github.com/merliseclyde/BAS/issues/40)

## Bug Fixes 

* Fixed valgrind error in src/ZS_approx_null_np.c for invalid write noted in CRAN checks

*Fixed function declarations identified with LTO

* Added `contrast=NULL` argument to `bas.lm` and `bas.glm` so that non-NULL contrasts do not
trigger warning in `model.matrix` as of R 3.6.0.  [Bug #44](https://github.com/merliseclyde/BAS/issues/44)

* Added check for sample size equal to zero due to subsetting or missing data
[Bug #37](https://github.com/merliseclyde/BAS/issues/37)

## Other 

* Put ORCID in quotes in author list (per R-dev changes)

