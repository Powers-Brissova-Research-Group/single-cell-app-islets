library(dplyr)
library(biomaRt)
library(readr)

# load data----
load("DATA/Islets2.Rda")

# extract Gene description
ref<-read_tsv("DATA/refdata-cellranger-GRCh38-1.2.0/Homo_sapiens.GRCh38.84_gene_annotation_table.txt")


# fetch gene id info from biomart
mart <- useEnsembl(biomart = "ensembl", 
                   host = "uswest.ensembl.org",
                   dataset = "hsapiens_gene_ensembl", 
                   version = "80")

genemap <- getBM( attributes = c("ensembl_gene_id",'hgnc_symbol', "description",'gene_biotype', 'chromosome_name', 'start_position', 'end_position'),
                  filters = "ensembl_gene_id",
                  values = ref$gene_id,
                  mart = mart,
                  useCache = FALSE)

#r <- setdiff(ref$gene_id, genemap$ensembl_gene_id) # what is in A that isn't in B, setdiff(A, B)       
#head(r)

# bind the seurat data back to original GRC38 genes, NOT biomart
idx <- match( rownames(Islets), ref$GeneSymbol)
Islets_ensembl <- ref[ idx, ]

# now, bind the biomart data
idx2 <- match( Islets_ensembl$gene_id, genemap$ensembl_gene_id) 
Islets_ensembl <- genemap[ idx2, ]


# and sort by gene symbol
Islets_ensembl <- Islets_ensembl %>% arrange(hgnc_symbol)

# fill in empty symbols with ensemblID (doing this does not let the app default to home tab)¯\_(ツ)_/¯
#Islets_ensembl$hgnc_symbol <- ifelse(Islets_ensembl$hgnc_symbol == "", Islets_ensembl$ensembl_gene_id, Islets_ensembl$hgnc_symbol)

# fill in empty symbols with "no symbol available"
Islets_ensembl$hgnc_symbol <- ifelse(Islets_ensembl$hgnc_symbol == "", "no symbol available", Islets_ensembl$hgnc_symbol)


#leave at least one empty hgnc (at least one empty hgnc_symbol helps the app default to home tab but too many empty hgnc_symbol outputs empty dropdown gene list to select)¯\_(ツ)_/¯
Islets_ensembl$hgnc_symbol<- ifelse(Islets_ensembl$ensembl_gene_id=="ENSG00000273338","",Islets_ensembl$hgnc_symbol)

# and save
write.csv(Islets_ensembl, "DATA/gene_annotation_sorted2.csv",row.names=FALSE)


#gene annotation that makes the app default to home tab but many genes listed does not output umap/vlnplots:

#Fetch gene id, gene description info from biomart
ensembl = useMart(
  "ensembl",
  host = "uswest.ensembl.org",
  dataset = "hsapiens_gene_ensembl" )
biomart_annotation <- getBM( attributes = c("ensembl_gene_id",'hgnc_symbol', "description",'gene_biotype', 'chromosome_name', 'start_position', 'end_position','source'),
                  filters = "hgnc_symbol",
                  values = rownames(Islets),
                  mart = ensembl,
                  useCache = FALSE)
write.csv(biomart_annotation, "DATA/biomart_annotation.csv")

