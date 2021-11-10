# BAS 1.6.0  Comments to CRAN

## submission to include  USE_FC_LEN_T for calls to Fortran subroutines with character strings in light of requirement under R 4.2.0


## Test environments

* local OS X install, R 4.1.2
* ubuntu 14.04 (on travis-ci), R 4.1.2 R-devel
* win-builder (r-devel)
* R-hub fedora-clang-devel (r-devel)
* R-hub solaris-x86-patched (r-patched)


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

