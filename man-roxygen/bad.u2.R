#' @param bad.u2 A numeric vector with elements belonging to \code{time}.  In
#' some cases, something goes wrong in the stratum, and the number of unmarked
#' fish captured should be ignored.  The values of \code{u2} in the entire row
#' will be set to NA for these strata. DO NOT SET the value of u2 to 0 because
#' this indicates that the trap was operating and captured no fish.
