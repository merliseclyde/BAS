# BAS 1.5.6  Comments to CRAN

resubmission to fix  errors found in package version 1.5.5: 
 
    CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2020-01-20 as check issues were still not
      corrected.
  
Fedora:   Internet access violation

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
* local fedora/clang fortran-9/valgrind R-devel (4.0.0)
* ubuntu 14.04 (on travis-ci), R 3.6.2 R-devel
* win-builder (r-devel)
* R-hub fedora-clang-devel (r-devel)
* R-hub solaris-x86-patched (r-patched)

Passes all checks on OSX/current, fedora/clang/R-devel, ubuntu/current/R-devel
Notes (windows) and error (Solaris) are false positives (see details below from checks).  Now issues identified via valgrind under fedora/clang fortran-9/valgrind R-devel (4.0.0)

In particular, the issues identified in the last submission on Fedora and Solaris have been identified and fixed for this submission

## R CMD check results for this submission

* On fedora-clang-devel (r-devel):  0 error | 0 warnings | 1 notes  

NOTE Maintainer: ‘Merlise Clyde <clyde@duke.edu>’
  
  
* On solaris-x86-patched (r-patched)  1 error | 0 warnings | 1 notes
  checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Merlise Clyde <clyde@duke.edu>’
  
    Found the following (possibly) invalid URLs:
    URL: https://dx.doi.org/10.1093/biomet/ass040
      From: man/bas.lm.Rd
            inst/doc/BAS-vignette.html
      Status: Error
      Message: libcurl error code 56:
        	OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 131

NOTE to CRAN:  url is valid,  but encounters time-out during check (passes check on all other platforms, and no issues with other urls in package)

* win-builder (r-devel)  0 errors | 0 warnings  | 1 notes 

  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Merlise Clyde <clyde@duke.edu>'
  
  New submission
  
  Package was archived on CRAN

 Possibly mis-spelled words in DESCRIPTION:
 
  Ghosh (40:12)
  Liang (25:30)
  Littman (40:22)
  Siow (24:34)
  Zellner (24:26)
  Zellner's (23:10)
  al (25:39)
  et (25:36)

NOTE to CRAN:  Notes are regarding maintainer (standard) and spelling (all correct)

## Reverse Dependencies

 
None

