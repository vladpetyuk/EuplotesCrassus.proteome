## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(message=FALSE)
# knitr::opts_chunk$set(echo=T, message=F, warning=F, fig.align='center', out.width='10cm')

## ----dataset_summary, echo=FALSE-----------------------------------------
library(knitr)
dat_sum <- data.frame(
    `Dataset Prefix`=c('Euplotes_1_SCX','Euplotes_1_HPRP_1','Euplotes_1_HPRP_2'),
    `Digestion Enzyme`=c('trypsin','trypsin','Glu-C (pH 7.5)'),
    `Fractionation Chromatrography Type`=c('SCX','HPRP','HPRP'),
    stringsAsFactors = FALSE,
    check.names = FALSE)
kable(dat_sum)

## ------------------------------------------------------------------------
fpath <- system.file("extdata",
                     "MSGFDB_GluC_StatCysAlk_10ppmParTol.txt",
                     package="EuplotesCrassus.proteome")
cat(readLines(fpath, n=12), sep = '\n')

## ----fetch_datasets, message=FALSE---------------------------------------
library(rpx)
id <- "PXD004333"
px <- PXDataset(id)
repoFiles <- pxfiles(px)
mzids <- grep('*msgfplus.mzid.gz', repoFiles, value=T)
system.time(pxget(px, mzids)) 

## ---- message=FALSE------------------------------------------------------
library(Biostrings)
fasta_clean <- readAAStringSet(
    system.file("extdata", 
                "Euplotes_Crassus_frameshifts.fasta",
                package="EuplotesCrassus.proteome"),
    format="fasta", nrec=-1L, skip=0L, use.names=TRUE)
fasta_marks <- readAAStringSet(
    system.file("extdata", 
                "Euplotes_Crassus_frameshifts_with_mark.fasta",
                package="EuplotesCrassus.proteome"),
    format="fasta", nrec=-1L, skip=0L, use.names=TRUE)
length(fasta_clean)

## ------------------------------------------------------------------------
library(MSnID)
trypscx <- grep('Euplotes_1_SCX_.*msgfplus.mzid.gz', repoFiles, value=T)
trypscxPrj <- MSnID()
system.time(trypscxPrj <- read_mzIDs(trypscxPrj, trypscx, backend = 'mzR'))

## ------------------------------------------------------------------------
trypscxPrj <- assess_termini(trypscxPrj, validCleavagePattern="[KR]\\.[^P]")
trypscxPrj <- apply_filter(trypscxPrj, "numIrregCleavages == 0")

## ------------------------------------------------------------------------
#' Rule on how to split the names.
#' Contig + Contaminants - main piece
#' comp - sequences with frameshifts
trypscxPrj.main <- apply_filter(trypscxPrj, "!grepl('comp', accession)")
trypscxPrj.fmsh <- apply_filter(trypscxPrj, "grepl('comp', accession)")
#' if peptide matches to the main piece we don't care about it
trypscxPrj.fmsh <- apply_filter(trypscxPrj.fmsh, 
                                "!(peptide %in% peptides(trypscxPrj.main))")
show(trypscxPrj.fmsh)

## ------------------------------------------------------------------------
trypscxPrj.fmsh$mme.ppm <- abs(mass_measurement_error(trypscxPrj.fmsh))
trypscxPrj.fmsh$score <- -log10(trypscxPrj.fmsh$`MS.GF.SpecEValue`)
trypscxPrj.fmsh <- apply_filter(trypscxPrj.fmsh, "mme.ppm < 10")

filtr <- MSnIDFilter(trypscxPrj.fmsh)
filtr$mme.ppm <- list(comparison="<", threshold=5.0)
filtr$score <- list(comparison=">", threshold=8.0)
#' pre-optimization with brute-force approach
filtr.grid <- optimize_filter(filtr, trypscxPrj.fmsh, fdr.max=0.05,
                              method="Grid", level="peptide", n.iter=20000)
evaluate_filter(trypscxPrj.fmsh, filtr.grid)
#' fine tune with optimization using simulated annealing technique
filtr.sann <- optimize_filter(filtr.grid, trypscxPrj.fmsh, fdr.max=0.05,
                              method="SANN", level="peptide", n.iter=20000)
evaluate_filter(trypscxPrj.fmsh, filtr.sann)
trypscxPrj.fmsh <- apply_filter(trypscxPrj.fmsh, filtr.sann)
show(trypscxPrj.fmsh)

## ---- message=FALSE------------------------------------------------------
#' extract only those that map frameshift sites
library(dplyr)
pepSeq <- unique(trypscxPrj.fmsh$pepSeq)
pepSeqMapped_to_clean <- pepSeq %>% 
    sapply(grep, x=fasta_clean) %>% 
    sapply(length) %>% 
    subset(.>0) %>% 
    names
pepSeqMapped_to_with_marks <- pepSeq %>% 
    sapply(grep, x=fasta_marks) %>% 
    sapply(length) %>% 
    subset(.>0) %>% 
    names
pepSeqFmsh_trypscx <- setdiff(pepSeqMapped_to_clean, pepSeqMapped_to_with_marks)
print(pepSeqFmsh_trypscx)

## ---- results='asis'-----------------------------------------------------
meta_tryp_scx <- trypscxPrj.fmsh %>%
    apply_filter('pepSeq %in% pepSeqFmsh_trypscx') %>%
    psms %>%
    select(spectrumFile,MS.GF.SpecEValue,mme.ppm,spectrumID,chargeState,peptide) %>%
    rename(SpecEValue = MS.GF.SpecEValue, charge = chargeState, `MME (ppm)`=mme.ppm) %>%
    mutate(spectrumFile = sub('_msgfplus.mzid.gz','',spectrumFile))
library(xtable)
print(xtable(meta_tryp_scx, display = c('d','s','e','f','s','d','s')),
      include.rownames=FALSE,
      comment = FALSE,
      size='scriptsize',
      floating = F)

## ------------------------------------------------------------------------
library(MSnID)
tryphprp <- grep('Euplotes_1_HPRP_1_.*msgfplus.mzid.gz', repoFiles, value=T)
tryphprpPrj <- MSnID()
system.time(tryphprpPrj <- read_mzIDs(tryphprpPrj, tryphprp, backend = 'mzR'))

## ------------------------------------------------------------------------
tryphprpPrj <- assess_termini(tryphprpPrj, validCleavagePattern="[KR]\\.[^P]")
tryphprpPrj <- apply_filter(tryphprpPrj, "numIrregCleavages == 0")

## ------------------------------------------------------------------------
tryphprpPrj.main <- apply_filter(tryphprpPrj, "!grepl('comp', accession)")
tryphprpPrj.fmsh <- apply_filter(tryphprpPrj, "grepl('comp', accession)")
tryphprpPrj.fmsh <- apply_filter(tryphprpPrj.fmsh, 
                                "!(peptide %in% peptides(tryphprpPrj.main))")
show(tryphprpPrj.fmsh)

## ------------------------------------------------------------------------
tryphprpPrj.fmsh$mme.ppm <- abs(mass_measurement_error(tryphprpPrj.fmsh))
tryphprpPrj.fmsh$score <- -log10(tryphprpPrj.fmsh$`MS.GF.SpecEValue`)
tryphprpPrj.fmsh <- apply_filter(tryphprpPrj.fmsh, "mme.ppm < 10")

filtr <- MSnIDFilter(tryphprpPrj.fmsh)
filtr$mme.ppm <- list(comparison="<", threshold=5.0)
filtr$score <- list(comparison=">", threshold=8.0)
filtr.grid <- optimize_filter(filtr, tryphprpPrj.fmsh, fdr.max=0.05,
                              method="Grid", level="peptide", n.iter=20000)
evaluate_filter(tryphprpPrj.fmsh, filtr.grid)
filtr.sann <- optimize_filter(filtr.grid, tryphprpPrj.fmsh, fdr.max=0.05,
                              method="SANN", level="peptide", n.iter=20000)
evaluate_filter(tryphprpPrj.fmsh, filtr.sann)
tryphprpPrj.fmsh <- apply_filter(tryphprpPrj.fmsh, filtr.sann)
show(tryphprpPrj.fmsh)

## ---- message=FALSE------------------------------------------------------
library(dplyr)
pepSeq <- unique(tryphprpPrj.fmsh$pepSeq)
pepSeqMapped_to_clean <- pepSeq %>% 
    sapply(grep, x=fasta_clean) %>% 
    sapply(length) %>% 
    subset(.>0) %>% 
    names
pepSeqMapped_to_with_marks <- pepSeq %>% 
    sapply(grep, x=fasta_marks) %>% 
    sapply(length) %>% 
    subset(.>0) %>% 
    names
pepSeqFmsh_tryphprp <- setdiff(pepSeqMapped_to_clean, pepSeqMapped_to_with_marks)
print(pepSeqFmsh_tryphprp)

## ---- results='asis'-----------------------------------------------------
meta_tryp_hprp <- tryphprpPrj.fmsh %>%
    apply_filter('pepSeq %in% pepSeqFmsh_tryphprp') %>%
    psms %>%
    select(spectrumFile,MS.GF.SpecEValue,mme.ppm,spectrumID,chargeState,peptide) %>%
    rename(SpecEValue = MS.GF.SpecEValue, charge = chargeState, `MME (ppm)`=mme.ppm) %>%
    mutate(spectrumFile = sub('_msgfplus.mzid.gz','',spectrumFile))
library(xtable)
print(xtable(meta_tryp_hprp, display = c('d','s','e','f','s','d','s')),
      include.rownames=FALSE,
      comment = FALSE,
      size='scriptsize',
      floating = F)

## ------------------------------------------------------------------------
library(MSnID)
gluchprp <- grep('Euplotes_1_HPRP_2_.*msgfplus.mzid.gz', repoFiles, value=T)
gluchprpPrj <- MSnID()
system.time(gluchprpPrj <- read_mzIDs(gluchprpPrj, gluchprp, backend = 'mzR'))

## ------------------------------------------------------------------------
gluchprpPrj <- assess_termini(gluchprpPrj, validCleavagePattern="E\\.[^P]")
gluchprpPrj <- apply_filter(gluchprpPrj, "numIrregCleavages == 0")

## ------------------------------------------------------------------------
gluchprpPrj.main <- apply_filter(gluchprpPrj, "!grepl('comp', accession)")
gluchprpPrj.fmsh <- apply_filter(gluchprpPrj, "grepl('comp', accession)")
gluchprpPrj.fmsh <- apply_filter(gluchprpPrj.fmsh, 
                                "!(peptide %in% peptides(gluchprpPrj.main))")
show(gluchprpPrj.fmsh)

## ------------------------------------------------------------------------
gluchprpPrj.fmsh$mme.ppm <- abs(mass_measurement_error(gluchprpPrj.fmsh))
gluchprpPrj.fmsh$score <- -log10(gluchprpPrj.fmsh$`MS.GF.SpecEValue`)
gluchprpPrj.fmsh <- apply_filter(gluchprpPrj.fmsh, "mme.ppm < 10")

filtr <- MSnIDFilter(gluchprpPrj.fmsh)
filtr$mme.ppm <- list(comparison="<", threshold=5.0)
filtr$score <- list(comparison=">", threshold=8.0)
filtr.grid <- optimize_filter(filtr, gluchprpPrj.fmsh, fdr.max=0.05,
                              method="Grid", level="peptide", n.iter=20000)
evaluate_filter(gluchprpPrj.fmsh, filtr.grid)
filtr.sann <- optimize_filter(filtr.grid, gluchprpPrj.fmsh, fdr.max=0.05,
                              method="SANN", level="peptide", n.iter=20000)
evaluate_filter(gluchprpPrj.fmsh, filtr.sann)
gluchprpPrj.fmsh <- apply_filter(gluchprpPrj.fmsh, filtr.sann)
show(gluchprpPrj.fmsh)

## ---- message=FALSE------------------------------------------------------
library(dplyr)
pepSeq <- unique(gluchprpPrj.fmsh$pepSeq)
pepSeqMapped_to_clean <- pepSeq %>% 
    sapply(grep, x=fasta_clean) %>% 
    sapply(length) %>% 
    subset(.>0) %>% 
    names
pepSeqMapped_to_with_marks <- pepSeq %>% 
    sapply(grep, x=fasta_marks) %>% 
    sapply(length) %>% 
    subset(.>0) %>% 
    names
pepSeqFmsh_gluchprp <- setdiff(pepSeqMapped_to_clean, pepSeqMapped_to_with_marks)
print(pepSeqFmsh_gluchprp)

## ---- results='asis'-----------------------------------------------------
meta_gluc_hprp <- gluchprpPrj.fmsh %>%
    apply_filter('pepSeq %in% pepSeqFmsh_gluchprp') %>%
    psms %>%
    select(spectrumFile,MS.GF.SpecEValue,mme.ppm,spectrumID,chargeState,peptide) %>%
    rename(SpecEValue = MS.GF.SpecEValue, charge = chargeState, `MME (ppm)`=mme.ppm) %>%
    mutate(spectrumFile = sub('_msgfplus.mzid.gz','',spectrumFile))
library(xtable)
print(xtable(meta_gluc_hprp, display = c('d','s','e','f','s','d','s')),
      include.rownames=FALSE,
      comment = FALSE,
      size='scriptsize',
      floating = F)

## ---- echo=FALSE, results='asis'-----------------------------------------
meta_tryp_scx$experiment <- "trypsin/SCX"
meta_tryp_hprp$experiment <- "trypsin/HPRP"
meta_gluc_hprp$experiment <- "Glu-C/HPRP"
meta <- rbind(meta_tryp_scx, meta_tryp_hprp, meta_gluc_hprp)
print(xtable(meta, display = c('d','s','e','f','s','d','s','s')),
      include.rownames=FALSE,
      comment = FALSE,
      size='tiny',
      floating = F)

## ---- echo=FALSE, results='asis'-----------------------------------------
toLatex(sessionInfo())

