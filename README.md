Data analysis part of the pipeline for confirming 
frameshifted protein sequences in *Euplotes crassus* with 
bottom-up LC-MS/MS proteomics.

To run install R (and preferrably RStudio).
```r
# install devtools if necessary
# note, devtools comes in a packages with RStudio
install.packages("devtools")
library("devtools")

install_github("vladpetyuk\EuplotesCrassus.proteome")
library("EuplotesCrassus.proteome")
vignette("euplotes_frameshifts")
```
