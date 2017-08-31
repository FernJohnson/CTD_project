#Using biomaRt to retrieve gene symbols if they are not given in annotation

library(GEOquery)
library(biomaRt)

GSE81071_eSet<-getGEO("GSE81071", GSEMatrix = TRUE)
GSE81071_eSet<-GSE81071_eSet[[1]]

#the platform has the entrez gene id for probes, but not the gene symbol

platf<-getGEO(annotation(GSE81071_eSet), AnnotGPL = TRUE)
UG<-attr(dataTable(platf), "table")[, "ENTREZ_GENE_ID"]
head(UG)

ensembl<- useMart("ensembl",dataset="hsapiens_gene_ensembl")

UGGS<-getBM(attributes = c('entrezgene', 'hgnc_symbol'),
            filters = 'entrezgene', 
            values = UG, 
            mart = ensembl)

head(UGGS)