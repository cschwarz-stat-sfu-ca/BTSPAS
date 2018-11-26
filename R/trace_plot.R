#' Creates trace plots of specified parameters showing the multiple chains and
#' the value of Rhat
#' 
#' Takes the MCMC object returned from a split and produces trace_plots for the
#' listed parameters. It shows a separate line on the plot for each chain and
#' also shows the value of Rhat
#' 
#' 
#' @param title A character string used for a title on reports and graphs
#' @param results The MCMC object containing the results from the call to
#' WinBugs/OpenBugs
#' @param parms_to_plot A character vector of names of parameters to plot.
#' These must match exactly to the parameter names used in the simulation.
#' @param panels How many plots to put on a page. It used split.screen to
#' format the plots.
#' @param mai Margins (inches) from the 4 borders. See par().
#' @param cex Character expansion factor. See par() for more details.
#' @return Nothing returned. Creates the plots
#' @author Bonner, S.J. \email{s.bonner@@stat.ubc.ca} and Schwarz, C. J.
#' \email{cschwarz@@stat.sfu.ca}
#' @keywords ~models ~trace plot
#' @examples
#'  
#' \dontrun{
#' # Create trace plots of the logitP parameters
#' # 
#' # extract the names of the variables from the MCMC object
#' varnames <- names(results$sims.array[1,1,]) 
#' # get the parms that start with logitP
#' parm.names <- varnames[grep("^logitP", varnames)] 
#' # create a pdf file of the plots
#' pdf(file="trace-logitP.pdf",sep=""))
#' trace_plot(title=title, results=results, 
#'     parms_to_plot=parm.names, panels=c(3,2))
#' dev.off()
#' 
#' # Or if you want an interactive display
#' par(ask=TRUE)
#' trace_plot(title=title, results=results, 
#'     parms_to_plot=parm.names, panels=c(2,1))
#' } % end of dontrun
#' 
trace_plot <- function(title=" ", results=NULL, parms_to_plot=NULL, panels=c(1,1),
   mai=if(prod(panels)>1){c(.4,.4,.4,.4)} else {c(1.02,0.82,0.82,0.42)},
   cex=if(prod(panels)>1){.5} else {1}
   ) {
#
# Takes the MCMC object from the fit (could be TPSDE etc), a list of parameters and produces
# the traceplots.
#   
# title - title of the plot
# results - the MCMC object containing the necessary information
# parms_to_plot - character vector containing the names of the parms to plot
#                e.g. c("logitP[1]", "logitP[2]")
#               - this should be an exact match
# panels  - how the plotting page should be arranged
# mai      - how big to make the margins around the plots
# cex      - character expansion factor
#

varnames <- names(results$sims.array[1,1,])

index <- match(parms_to_plot, varnames) # find where these parms exist in the array
plots_per_page <- prod(panels)

for(i in seq(1,length(index),plots_per_page)){  # plot potentially multiple plots/page
   split.screen(panels)
   par(new=TRUE)
   for(j in i:min(i+plots_per_page-1,length(index))){
      screen(j-i+1)
      par(mai=mai)
      par(cex=cex)
      matplot(results$sims.array[,,index[j]], type="l", 
              main=paste(title), lty=1,
              ylab='Estimate')
      text(x=1, y=max(results$sims.array[,,index[j]]), 
          label=varnames[index[j]], adj=c(0,1), cex=if(prod(panels)>1){2} else {1})
      if(results$n.chains >1){   # Only print Rhat if number of chains is >1
         text(x=dim(results$sims.array)[1]  ,y=max(results$sims.array[,,index[j]]), 
          label=paste("Rhat=",round(results$summary[index[j],"Rhat"],1)), adj=c(1,1), 
          cex=if(prod(panels)>1){2} else {1})
      }
      if(results$n.chains == 1){
         text(x=dim(results$sims.array)[1]  ,y=max(results$sims.array[,,index[j]]), 
          label="No RHat avail", adj=c(1,1), 
          cex=if(prod(panels)>1){2} else {1})
      }
   }
   close.screen(all.screens=TRUE)
}

return(NULL) # nothing to return
} # end of function
