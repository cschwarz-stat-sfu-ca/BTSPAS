#' @param logitP.fixed A numeric vector (could be null) of the time strata
#' where the logit(P) would be fixed. Typically, this is used when the capture
#' rates for some strata are 0 and logit(P) is set to -10 for these strata. The
#' fixed values are given in \code{logitP.fixed.values}
#' @param logitP.fixed.values A numerical vector (could be null) of the fixed
#' values for logit(P) at strata given by logitP.fixed. Typically this is used
#' when certain strata have a 0 capture rate and the fixed value is set to -10
#' which on the logit scale gives p[i] essentially 0. Don't specify values such
#' as -50 because numerical problems could occur in JAGS.
