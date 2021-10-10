
## Major change

* Fixed issue with reverse dependency on actuar package with sd() and var() functions as
notified by R team.


## Test environments
* local OS X install, R 4.1.1
* local Windoze 10, R 4.1.1

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 1 NOTES:

* checking installed package size ... NOTE
  installed size is  5.3Mb
  sub-directories of 1Mb or more:
    doc   4.9Mb

* Found the following (possibly) invalid URLs:
  URL: https://www.jstor.org/stable/2332748
    From: inst/doc/e-Bias-from-incomplete-sampling.html
    Status: 403
    Message: Forbidden
    
    ***CJS*** URL is valid.
    
## Reverse dependencies
None.