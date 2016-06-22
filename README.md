Data analysis part of the pipeline for confirming 
frameshifted protein sequences in *Euplotes crassus* with 
bottom-up LC-MS/MS proteomics. This is supplementary
Bioconductor/R packages for the manuscript:

> **"Widespread Abrogation of Triplet Translation Continuity and 
> Stop Codon Function in Euplotes"**

> Alexei V. Lobanov, Stephen M. Heaphy, Anton A. Turanov, 
> Maxim V. Gerashchenko, Sandra Pucciarelli, Raghul R. Devaraj, 
> Fang Xie, Vladislav A. Petyuk, Richard D. Smith, 
> Lawrence A. Klobutcher, John F. Atkins, Cristina Miceli, Dolph L. Hatfield, 
> Pavel V. Baranov, Vadim N. Gladyshev


To run install R (and preferrably RStudio).
```r
# install devtools if necessary
# note, devtools comes in a packages with RStudio
install.packages("devtools")
library("devtools")

# add path to Bioconductor repositories
source("http://bioconductor.org/biocLite.R")
options(repos=biocinstallRepos(character()))

install_github("vladpetyuk/EuplotesCrassus.proteome", build_vignettes=TRUE)

library("EuplotesCrassus.proteome")
vignette("euplotes_frameshifts")
```
