# Security Policy

## Supported Versions

Supported with security updates.

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Reporting a Vulnerability

Please submit any vulnerability reports under [Github Issues](https://github.com/merliseclyde/bark/issues) and maintainers will address as soon as possibl

## Expectations

This package utilizes C code for efficiency and allocates/frees memory. 
The package is checked for memory leaks prior to releases to CRAN using 
ASAN/UBSBAN. The package is distributed via CRAN   https://CRAN.R-project.org/package=bark which reports additional checks. The development version may be installed from GitHub https://github.com/merliseclyde/bark which is checked via github actions 
(users may check that the current version has a passing badge before installing)
Bugs are reported via the Issue tracker and handled as soon as possible.  
(See link above)

## Assurance

It is highly unlikely that malicious code would be added to the package. All submissions to CRAN require verification via the maintainer's email, which is protected via two factor authentication.  Any pull requests for contributions on github are verified by the lead maintainer.  Based on the Code of Conduct and Contributing Guidelines any modifications should include unit tests that cover the additional code blocks.