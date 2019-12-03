# render all my vignettes into html in the vignette directory

library(plyr)
library(rmarkdown)
# navigate to the vignettes directory
files <- dir()
files <- files[ grepl("Rmd$", files)]

plyr::l_ply(files, function(x){
  cat("Starting to render ", x, "\n")
    rmarkdown::render(x)
})