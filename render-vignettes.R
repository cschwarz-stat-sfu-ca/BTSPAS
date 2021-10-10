# render all my vignettes into html in the vignette directory
# navigate to the vignettes directory first, before running this job

run.parallel <- FALSE

library(plyr)
library(rmarkdown)
if(run.parallel){
  library(doMC)

  ## Specify the number of cores
  registerDoMC(6)

  ## Check how many cores we are using
  getDoParWorkers()
}


files <- dir()
files <- files[ grepl("Rmd$", files)]
files

plyr::l_ply(files[5], function(x){
  cat("Starting to render ", x, "\n")
    rmarkdown::render(x)
},  .parallel=run.parallel)

#rmarkdown::render(files[6])
