
## Major change

* None. Bug fix to deal with issue in Petersen package which was traced to 
an bug in the R2jags::jags() function when pD=TRUE but pD=0, it is not included
in the BUGSoutput object.

## Test environments
* local OS X install, R 4.4.1
* devtools::check_win_release()
* devtools::check_win_devel()

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 1 NOTES:

❯ checking installed package size ... NOTE
    installed size is  5.3Mb
    sub-directories of 1Mb or more:
      doc   4.9Mb

Some rhub checks failed because JAGS not installed in those machines.

## Reverse dependencies
### revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
 

