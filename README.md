#  Single cell gene expression atlas of human pancreatic islets 

This repository contains R scripts for R shiny web application to browse single cell gene expression data from Shrestha et al., bioRXiv(2021).This expansive dataset contains ~45,000 cells from 5 donor pancreatic islets and provides a comprehensive map of human islet biology. It will be an extremely useful resource that complements existing single cell datasets from human islets. 

## Installation

The web app can be accessed at https://powersbrissovalab.shinyapps.io/scRNAseq-Islets/ or the web app can run locally through app.R

## Note about data

The data (stored as `.rda` file) was too large for Github's single file limit, so it was chunked into smaller files using the 7zip utility. To re-generate the data file, please use 7zip to extract the data `DATA/Islets.Rda.zip.001` through `DATA/Islets.Rda.zip.004`. 7zip will extract to to `Islets.Rda`, which is what the R Shiny app expects. You may then delete the chunked files.

---
Page last edited on 2020-12-04

