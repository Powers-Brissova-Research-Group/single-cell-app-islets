
# load data----
load("DATA/Islets2.Rda")

# extract Gene description
ref<-read_tsv("DATA/refdata-cellranger-GRCh38-1.2.0/Homo_sapiens.GRCh38.84_gene_annotation_table.txt")


# fetch gene id info from biomart
library(biomaRt)
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

# fill in empty symbols with ensemblID
Islets_ensembl$hgnc_symbol <- ifelse(Islets_ensembl$hgnc_symbol == "", Islets_ensembl$ensembl_gene_id, Islets_ensembl$hgnc_symbol)

# and sort by gene symbol
Islets_ensembl <- Islets_ensembl %>% arrange(hgnc_symbol)

# and save... as you'd like
write.csv(Islets_ensembl, "DATA/gene_annotation_sorted.csv",row.names=FALSE)



