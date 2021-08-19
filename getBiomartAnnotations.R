
#load data----
load("DATA/Islets2.Rda")

#Extract Gene description
ref<-read_tsv("DATA/refdata-cellranger-GRCh38-1.2.0/Homo_sapiens.GRCh38.84_gene_annotation_table.txt")


#fetch gene id info from biomart
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", 
                   host = "uswest.ensembl.org",
                   dataset = "hsapiens_gene_ensembl", 
                   version = "80")

genemap <- getBM( attributes = c("ensembl_gene_id",'hgnc_symbol', "description",'gene_biotype', 'chromosome_name', 'start_position', 'end_position','source'),
                  filters = "ensembl_gene_id",
                  values = ref$gene_id,
                  mart = mart,
                  useCache = FALSE)

r <- setdiff(ref$gene_id, genemap$ensembl_gene_id) # what is in A that isn't in B, setdiff(A, B)       
head(r)

# bind the seurat data back to original GRC38 genes, NOT biomart
idx <- match( rownames(Islets), ref$GeneSymbol)
Islets_ensembl <- ref[ idx, ]

idx2 <- match( Islets_ensembl$gene_id, genemap$ensembl_gene_id) 
Islets_ensembl <- genemap[ idx, ]

# Fill in empty symbols
Islets_ensembl$hgnc_symbol <- ifelse(Islets_ensembl$hgnc_symbol == "", "No symbol available", Islets_ensembl$hgnc_symbol)
