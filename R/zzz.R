

#' Message to display when package is loaded
#' 
#' @keywords internal
#'

.onAttach <- function(libname,pkgname){

  packageStartupMessage("***** BTSPAS: Bayesian Time Stratified Petersen Analysis System - Version 2026-03-01 (2026.3.1) ***** \n\n",
        "      Help available with  help(package='BTSPAS') \n",
        '      Several vignettes are available. See browseVignettes(package="BTSPAS") \n\n')
}







