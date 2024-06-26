% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TimeStratPetersenDiagError_fit.R
\name{TimeStratPetersenDiagError_fit}
\alias{TimeStratPetersenDiagError_fit}
\title{Wrapper (*_fit) to call the Time Stratified Petersen Estimator
with Diagonal Entries function.}
\usage{
TimeStratPetersenDiagError_fit(
  title = "TSDPE",
  prefix = "TSPDE-",
  time,
  n1,
  m2,
  u2,
  sampfrac = rep(1, length(u2)),
  jump.after = NULL,
  bad.n1 = c(),
  bad.m2 = c(),
  bad.u2 = c(),
  logitP.cov = as.matrix(rep(1, length(n1))),
  logitP.fixed = NULL,
  logitP.fixed.values = NULL,
  n.chains = 3,
  n.iter = 2e+05,
  n.burnin = 1e+05,
  n.sims = 2000,
  tauU.alpha = 1,
  tauU.beta = 0.05,
  taueU.alpha = 1,
  taueU.beta = 0.05,
  prior.beta.logitP.mean = c(logit(sum(m2, na.rm = TRUE)/sum(n1, na.rm = TRUE)), rep(0,
    ncol(as.matrix(logitP.cov)) - 1)),
  prior.beta.logitP.sd = c(stats::sd(logit((m2 + 0.5)/(n1 + 1)), na.rm = TRUE), rep(10,
    ncol(as.matrix(logitP.cov)) - 1)),
  tauP.alpha = 0.001,
  tauP.beta = 0.001,
  run.prob = seq(0, 1, 0.1),
  debug = FALSE,
  debug2 = FALSE,
  InitialSeed = ceiling(stats::runif(1, min = 0, max = 1e+06)),
  save.output.to.files = TRUE,
  trunc.logitP = 15,
  set.browser = FALSE
)
}
\arguments{
\item{title}{A character string used for a title on reports and graphs}

\item{prefix}{A character string used as the prefix for created files. All
created graph files are of the form prefix-xxxxx.pdf.}

\item{time}{A numeric vector of time used to label the strata. For example,
this could be julian week for data stratified at a weekly level.}

\item{n1}{A numeric vector of the number of marked fish released in each
time stratum.}

\item{m2}{A numeric vector of the number of marked fish from n1 that are
recaptured in each time stratum. All recaptures take place within the
stratum of release.}

\item{u2}{A numeric vector of the number of unmarked fish captured in each
stratum. These will be expanded by the capture efficiency to estimate the
population size in each stratum.}

\item{sampfrac}{\strong{Deprecated} because it really doesn't work as intended.
You must remove all references to sampfrac from your code. Contact \email{cschwarz.stat.sfu.ca@gmail.com}
for more information.}

\item{jump.after}{A numeric vector with elements belonging to \code{time}.
In some cases, the spline fitting the population numbers should be allowed
to jump.  For example, the population size could take a jump when hatchery
released fish suddenly arrive at the trap.  The jumps occur AFTER the strata
listed in this argument.}

\item{bad.n1}{A numeric vector with elements belonging to \code{time}.  In
some cases, something goes wrong in the stratum, and the number of marked
fish releases should be discarded.  The values of \code{n1} and \code{m2} will be set to
0 for these strata.}

\item{bad.m2}{A numeric vector with elements belonging to \code{time}.  In
some cases, something goes wrong in the stratum, and the number of recovered
marked fish should be ignored. For example, poor handling is suspected to
induce handling induced mortality in the marked fish and so only very few
are recovered.  The values of \code{n1} and \code{m2} will be set to 0
for these strata.}

\item{bad.u2}{A numeric vector with elements belonging to \code{time}.  In
some cases, something goes wrong in the stratum, and the number of unmarked
fish captured should be ignored.  The values of \code{u2} in the entire row
will be set to NA for these strata. DO NOT SET the value of u2 to 0 because
this indicates that the trap was operating and captured no fish.}

\item{logitP.cov}{A numeric matrix for covariates to fit the
logit(catchability). Default is a single intercept, i.e. all strata have the
same mean logit(catchability).}

\item{logitP.fixed}{A numeric vector (could be null) of the time strata
where the logit(P) would be fixed. Typically, this is used when the capture
rates for some strata are 0 and logit(P) is set to -10 for these strata. The
fixed values are given in \code{logitP.fixed.values}}

\item{logitP.fixed.values}{A numerical vector (could be null) of the fixed
values for logit(P) at strata given by logitP.fixed. Typically this is used
when certain strata have a 0 capture rate and the fixed value is set to -10
which on the logit scale gives p[i] essentially 0. Don't specify values such
as -50 because numerical problems could occur in JAGS.}

\item{n.chains}{Number of parallel MCMC chains to fit.}

\item{n.iter}{Total number of MCMC iterations in each chain.}

\item{n.burnin}{Number of burn-in iterations.}

\item{n.sims}{Number of simulated values to keeps for posterior
distribution.}

\item{tauU.alpha}{One of the parameters along with \code{tauU.beta} for the
prior for the variance of the random noise for the smoothing spline.}

\item{tauU.beta}{One of the parameters along with \code{tauU.alpha} for the
prior for the variance of the random noise for the smoothing spline.}

\item{taueU.alpha}{One of the parameters along with \code{taueU.beta} for
the prior for the variance of noise around the spline.}

\item{taueU.beta}{One of the parameters along with \code{taueU.alpha} for
the prior for the variance of noise around the spline.}

\item{prior.beta.logitP.mean}{Mean of the prior normal distribution for
logit(catchability) across strata}

\item{prior.beta.logitP.sd}{SD of the prior normal distribution for
logit(catchability) across strata}

\item{tauP.alpha}{One of the parameters for the prior for the variance in
logit(catchability) among strata}

\item{tauP.beta}{One of the parameters for the prior for the variance in
logit(catchability) among strata}

\item{run.prob}{Numeric vector indicating percentiles of run timing should
be computed.}

\item{debug}{Logical flag indicating if a debugging run should be made. In
the debugging run, the number of samples in the posterior is reduced
considerably for a quick turn around.}

\item{debug2}{Logical flag indicated if additional debugging information is
produced. Normally the functions will halt at \code{browser()} calls to
allow the user to peek into the internal variables. Not useful except to
package developers.}

\item{InitialSeed}{Numeric value used to initialize the random numbers used
in the MCMC iterations.}

\item{save.output.to.files}{Should the plots and text output be save to the files
in addition to being stored in the MCMC object?}

\item{trunc.logitP}{Truncate logit(P) between c(=trunc.logitP, trunc.logitP) when plotting the
logitP over time. Actual values of logit(P) are not affected.}

\item{set.browser}{Should the function enter browser model when called (useful for debugging)}
}
\value{
An MCMC object with samples from the posterior distribution. A
series of graphs and text file are also created in the working directory.
}
\description{
Takes the number of marked fish released, the number of recaptures, and the
number of unmarked fish and uses Bayesian methods to fit a fit a spline
through the population numbers and a hierarchical model for the trap
efficiencies over time.  The output is written to files and an MCMC object
is also created with samples from the posterior.
}
\details{
Normally, the wrapper (*_fit) function is called which then calls the
fitting routine.

Use the \code{\link{TimeStratPetersenNonDiagError_fit}} function for cases
where recaptures take place outside the stratum of release.
}
\references{
Bonner, S. J., & Schwarz, C. J. (2011). 
Smoothing population size estimates for Time-Stratified Mark-Recapture experiments Using Bayesian P-Splines. 
Biometrics, 67, 1498-1507. 
\doi{10.1111/j.1541-0420.2011.01599.x}

Schwarz, C. J., & Dempson, J. B. (1994).
Mark-recapture estimation of a salmon smolt population. 
Biometrics, 50, 98-108.

Schwarz, C.J., D. Pickard, K. Marine and S.J. Bonner. 2009. Juvenile Salmonid Outmigrant
Monitoring Evaluation, Phase II - December 2009. Final Technical Memorandum for the Trinity River
Restoration Program, Weaverville, CA. 155 pp. + appendices available at
\url{https://www.trrp.net/library/document/?id=369}
}
\author{
Bonner, S.J. \email{sbonner6@uwo.ca} and 
Schwarz, C. J. \email{cschwarz.stat.sfu.ca@gmail.com}.
}
\keyword{~models}
\keyword{~smooth}
