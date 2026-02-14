library(Seurat)


load("DATA/Islets3.Rda")


# Find features that are not expressed in at least 10 cells
features <- rownames(Islets)[Matrix::rowSums(Islets[["RNA"]]@data > 0) >= 10]


#x <- Islets@assays$RNA@data[0:100, 0:100]


obj <- DietSeurat(
  Islets,
  layers = c("data"),
  features = features,
  assays = c("RNA"),
  dimreducs = c("umap"),
  graphs = NULL
)

dim(Islets)
dim(obj)


saveRDS(obj, "DATA/Islets4-slimmed.Rds")
