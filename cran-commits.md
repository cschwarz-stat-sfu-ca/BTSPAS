## Test environments
* local OS X install, R 3.6.1
* local Windoze 10, R 3.6.1

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTES:
  
* checking CRAN incoming feasibility ... NOTE
   Maintainer: ‘Carl James Schwarz <cschwarz.stat.sfu.ca@gmail.com>’
   Version contains large components (2020.1)  
  
  This is the 2020 release of the package.

## Reverse dependencies
None.

## Other notes
Vignettes are computationally intensive. Used the rmarkdown_notangle to prevent .R files from being created (and run) in the doc directory during package checks.