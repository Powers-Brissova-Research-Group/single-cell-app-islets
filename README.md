#  Single cell gene expression atlas of human pancreatic islets

This repository contains R scripts for R shiny web application to browse single cell gene expression data from Shrestha et al., JCI Insight(2021) doi:10.1172/jci.insight.151621. This expansive dataset contains ~45,000 cells from 5 donor pancreatic islets and provides a comprehensive map of human islet biology. It will be an extremely useful resource that complements existing single cell datasets from human islets.

## Installation & deployment

The web app can be accessed at https://powersbrissovalab.shinyapps.io/scRNAseq-Islets/ or the web app can run locally through app.R.

For AWS/Container deployment, see DEPLOYMENT.md

## Note about data

The single cell object is stored  `.Rds` file, which is too large for Github's single file limit, so it was chunked the 7zip utility. To re-generate the data file, please use 7zip to extract the data `Islets4-slimmed.Rds.gz.*` (multipart). 7zip will extract to `Islets4-slimmed.Rds`, which is what the R Shiny app expects. You may then delete the chunked files.

## Releases

### 2026-02-13

- Testing deployment to AWS Lightsail via ECR (replaciing shinyapps.io)
- Violin plots render without individual cell points by default
- Per-gene data caching (FetchData and DotPlot computations)
- Elapsed time logging on all plots and data fetches for diagnostics
- Consistent two-column layout across all plot pages: plots on left, options sidebar on right
- Gene info table transposed and moved to sidebar on all plot tabs

- **Violin Plot**

  - "Show individual cells" checkbox (off by default) with cell size slider

  - "Show boxplot" checkbox (on by default)

  - Resizable plots (height and width sliders)

- **UMAP**

  - Square aspect ratio (coord_fixed)

  - Resizable via single plot size slider

  - Adjustable point size

  - Color scheme selector: Grey-Red (default), Viridis, Magma, Inferno, Plasma, Cividis

- **Dot Plot**

  - Resizable via single plot size slider (height at 0.8 ratio of width)

- **Expression Values**

  - Adjustable decimal places (0-10, default 2)

- **Experimental Summary**

  - Simplified to plain static table (removed DT pagination/search)

### 2021-07-22

- Initial release on shinyapps.io
- Original Seurat object with dittoPlot-based visualizations

---
Page last edited on 2026-02-14
