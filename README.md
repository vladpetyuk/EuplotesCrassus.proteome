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

To run install:

1. R language to statistical computing: [R](https://cloud.r-project.org/)
2. IDE for R: [RStudio](https://www.rstudio.com/products/rstudio/download/)
3. `pdflatex` for compiling `pdf` files form LaTeX code. 
    The basic LaTeX installation should be sufficient. 
    Although during the R vignette compilation
    it may require installation a few extra packages. The packages are 
    installed on the fly after accepting by clicking `OK` button. 
    To install LaTeX follow these links:
    1. Windows: [MikTex](http://miktex.org/download)
    2. Mac OS: [MacTex](https://tug.org/mactex/)
    3. Linus: On Linux `pdflatex` is likely to be present, otherwise use
       your package manager. e.g.
       
       `sudo apt-get install texlive-latex-base`

Copy and paste the following code into R console. However, 
if you do not intend to re-compile the vignette (and thus reproduce the
analysis) set the `build_vignettes` to `FALSE`.

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

To reproduce the vignette compilation install the package with 
`build_vignettes=TRUE` or directly open and recompile 
`euplotes_frameshifts.Rmd` from the packages source.


