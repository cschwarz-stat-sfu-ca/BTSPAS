#' Compute percentiles of the run timing distribution.
#' 
#' Take the posterior sample of U[1,...nstrata] and compute the percentiles of the run
#' timing. 
#' This uses the quantile() function from the "actuar" package which is designed to compute
#' quantiles of grouped data.
#' 
#' @param time List of sample weeks. It is assumed that there are no salmon prior to the first value
#            in time, and after the last value in time
#' @param U matrix of posterior samples. Each row is a sample from the posterior.
#          Columns correspond to U[1]...U[nstrata]
#' @param prob Quantiles of the run timing to estimate. 

#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the working directory.
#' @author Bonner, S.J. \email{s.bonner@@stat.ubc.ca} and Schwarz, C. J.
#' \email{cschwarz@@stat.sfu.ca}
#' @references Refer to the Trinity River Restoration Project report by
#' Schwarz, C.J. et al. (2009) available at
#' \url{http://www.stat.sfu.ca/~cschwarz/Consulting/Trinity/Phase2}. Please
#' contact \email{cschwarz@stat.sfu.ca} for more details. %% ~put references to
#' the literature/web site here ~
#
#' @export RunTime
#' @import plyr
#' @importFrom actuar grouped.data
# The actuar pacakge cm() function conflicts with another package.
# The following excludes is
# See https://stackoverflow.com/questions/51899220/import-all-the-functions-of-a-package-except-one-when-building-a-package
#' @rawNamespace import(actuar, except = cm) 

# 2018-12-14 CJS converted from a for() loop to adply()

RunTime <- function(time, U, prob=seq(0,1,.1)) {
  timing <- c(min(time):(1+max(time)))
  q.U <- plyr::adply(U, 1, function(U.sample, timing){
       quant <- quantile(grouped.data(Group=timing, Frequency=U.sample), prob=prob)
       quant
  }, timing=timing, .id=NULL)
  q.U
}
