# 2020-12-15 CJS Fixed problem with u2==NA in the plots
# 2015-06-10 CJS Updated to produce ggplot2 plots
# 2014-09-21 CJS change all Inf to NA
# 2012-01-22 CJS made X/Y axis limits the same so p-value prints properly
# 2011-06-13 CJS returned bayesian p-values

PredictivePosteriorPlot.TSPDE.WHCH <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the discrep matrix has the following columns
#     (1-2)  o,s Freeman-Tukey measures for m2
#     (3-4)  o,s Freeman-Tukey measures for u2.A
#     (5-6)  o,s Freeman-Tukey measures for u2.N
#     (7-8)  o,s Freeman-Tukey for m2+u2.A+u2.N

# Change any Inf to NA
temp <- is.infinite(discrep) & !is.na(discrep)
if(sum(temp, na.rm=TRUE)>0){cat(sum(temp, na.rm=TRUE), " infinite discrepancy measures set to NA\n")}
discrep[ temp ] <- NA

#browser()
titles <- c("Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.A", 
            "Freeman-Tukey for u2.N", 
            "Total Freeman-Tukey")

saved_p_values <- rep(NA, length(titles))

discrep.df <- data.frame(discrep)
plot.list<- plyr::llply(1:4, function(i){
  p.value <- mean(discrep[,2*i-1]<discrep[,2*i],na.rm=TRUE)
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

bigplot <- do.call(arrangeGrob, c(plot.list, list(ncol=2)))

gof <- list(bp.plot=bigplot,  bp.values=data.frame(test.names=titles, p.value=saved_p_values, stringsAsFactors=FALSE))
gof
}
