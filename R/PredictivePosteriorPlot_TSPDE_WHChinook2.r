#' @rdname PredictivePosterior.TSPDE


# 2015-06-10 CJS Converted to ggplot()
# 2014-09-01 CJS Change Inf to NA 
# 2012-01-22 CJS Made X/Y axis limits the same so that Bayesian p-value prints properly
# 2011-06-13 CJS returned p-values
# 2010-03-29 CJS First creation of routine

PredictivePosteriorPlot.TSPDE.WHCH2 <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the discrepancy measure has 16 columns for the Bayesian p-value plot
#     ( 1- 2)  o,s Freeman-Tukey measures for m2
#     ( 3- 4)  o,s Freeman-Tukey measures for u2.A.YoY
#     ( 5- 6)  o,s Freeman-Tukey measures for u2.N.YoY
#     ( 7- 8)  o,s Freeman-Tukey measures for u2.A.1
#     ( 9-10)  o,s Freeman-Tukey measures for u2.N.1
#     (11-12)  o,s Freeman-Tukey for u2.A.YoY+u2.N.YoY
#     (13-14)  o,s Freeman-Tukey for u2.A.1  +u2.N.1
#     (15-16)  o,s Freeman-Tukey for all data (m2, YoY and Age 1)`

# Change any Inf to NA
temp <- discrep == Inf | discrep == -Inf
if(sum(temp)>0){cat(sum(temp), " infinite discrepancy measures set to NA\n")}
discrep[ temp ] <- NA

#browser()
titles <- c("Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.A.YoY", 
            "Freeman-Tukey for u2.N.YoY", 
            "Freeman-Tukey for u2.A.1", 
            "Freeman-Tukey for u2.N.1", 
            "Freeman-Tukey for YoY",
            "Freeman-Tukey for Age 1",
            "Total Freeman-Tukey")
saved_p_values <- rep(NA, length(titles))

discrep.df <- data.frame(discrep)
plot.list<- llply(1:8, function(i){
  p.value <- sum(discrep[,2*i-1]<discrep[,2*i],na.rm=TRUE)/nrow(discrep)
  saved_p_values[i] <<- p.value
  
  bp.plot <- ggplot(data=discrep.df, aes_string(x=colnames(discrep.df)[2*i], y=colnames(discrep.df)[2*i-1]))+
    ggtitle(titles[i])+
    geom_point()+
    xlab("Simulated")+ylab("Observed")+
    geom_abline(intercept=0, slope=1)+
    annotate("text", x=Inf,y=-Inf, hjust=1, vjust=0,
               label=paste("Bayesian GOF P:",formatC(p.value, digits=2, format="f")))
  bp.plot 
})

bigplot <- do.call(marrangeGrob, c(plot.list, list(ncol=2, nrow=2, top=NULL)))
gof <- list(bp.plot=bigplot,  bp.values=data.frame(test.names=titles, p.value=saved_p_values, stringsAsFactors=FALSE))

gof
}
