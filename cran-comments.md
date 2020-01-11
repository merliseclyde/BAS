#  BAS 1.5.4 Comments to CRAN  (resubmission)
 
This submission addresses additional issues found on the *last* version released on CRAN:
 LTO <https://urldefense.proofpoint.com/v2/url?u=https-3A__www.stats.ox.ac.uk_pub_bdr_LTO_BAS.out&d=DwICAw&c=imBPVzF25OnBgGmVOlcsiEgHoG1i6YHLR0Sj_gZ4adc&r=NOkxkvdFOOffXzeTY2kgZQ&m=52W89JkPbOpG86B0ILRzgnU71_hsgNFiQNKrlb5-Nt0&s=r-lKhBavxLKaoTr_9AlOwdBxfkFOR0sTSi-oR-WKzxU&e= >
 valgrind <https://urldefense.proofpoint.com/v2/url?u=https-3A__www.stats.ox.ac.uk_pub_bdr_memtests_valgrind_BAS&d=DwICAw&c=imBPVzF25OnBgGmVOlcsiEgHoG1i6YHLR0Sj_gZ4adc&r=NOkxkvdFOOffXzeTY2kgZQ&m=52W89JkPbOpG86B0ILRzgnU71_hsgNFiQNKrlb5-Nt0&s=xZMf4TAeqaSYWV8VES3z6OGZLkp320XviDiKS8Wj8Yg&e= >

Subsequent testing on fedora with clang (locally) and debian/gcc/valgrind via rhub indicate that the package passes above checks for issues identified above.  

The ERRORS in the *last* submission for Solaris remain (errors were allowed for last submission).

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

