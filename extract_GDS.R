#Extracting data from GDS (GEO datasets) accessions


library(GEOquery)

#First identify GDS accessions

geogds<-c()
for(i in geoaccs)
{
  
  
  if(grepl("GDS", i, fixed = TRUE))
    geogds<-c(geogds,i)
  
}


for(i in geogds)
{
#Download GDS - replace i with string of acession.

  gds<-getGEO(i)


#Convert to eSet - data is not log transformed.
  eset<-GDS2eSet(gds, do.log2 = FALSE)

#Assign name to eSet matching accesion.

  assign(paste(i, "eset",sep="_"),eset)
  
}

#Now you have all the expression sets, the expression data can be extracted
#and annotated with gene symbols

#In most cases the esets annotation slot is empty - but it can be found on
#the GEO website

#Download platform information, extract probe IDs and corresponding gene symbols
platf<-getGEO("GPL570", AnnotGPL = TRUE)
IDs<-attr(dataTable(platf), "table")[, c("ID", "Gene symbol")]

#Extract expression values and probe IDs (rownames)
exGDS3940<-exprs(GDS3940_eset)
rows<-rownames(exGDS3940)

#Find the probe IDs that correspond to a gene symbol

symb<-rep(0, length(rows))
for(i in 1:length(rows))
{
  symb[i]<-as.character(IDs[which(IDs[,1] == rows[i]), 2])
}

rownames(exGDS3940)<-symb


#Remove rows without a gene symbol.

sel<-apply(as.matrix(symb, ncol=1), 1, function(x){if(x=="") return(FALSE) else return(TRUE)})

selExpr<-c()
selSymb<-c()

for(i in 1:nrow(exGDS3940))
{
  if(sel[i])
  {
    selExpr<-rbind(selExpr, exGDS3940[i,])
    selSymb<-rbind(selSymb, symb[i])
  }
}

rownames(selExpr)<-selSymb

#Aggregate repeated measurements of the same gene symbol by their mean

selExprm<-aggregate(selExpr, by=list(row.names(selExpr)), FUN=mean)

#Aggregating sorts rows alphabetically and removes their rownames
#So adjust rownames to match, and assign

genenames<-unique(rownames(selExpr))

genenames<-sort(genenames)

rownames(selExprm)<-genenames

#Name the processed data with it's accesion number

assign(paste("GDS3940", "exp", sep = "_"), selExprm)

#Remove the first column - introduced by aggregate, no info
#Write table to CSV

GDS3940_exp<-GDS3940_exp[,-1]

write.csv(GDS5421_exp, file="GDS5421_exp.csv")

#Transform by ranks

GDS3940_exprank<-apply(GDS3940_exp, 2, rank)

write.csv(GDS3940_exprank, file="GDS3940_exprank.csv")

#Extract pheno data from eset
#Transform so samples accessions are column headings, to match exp data

phen<-t(pData(phenoData(GDS3940_eset)))

write.csv(phen, file="GDS3940_pheno.csv")


#Alternatively the ExtractExpX function can be used, rather than running all this code multiple times
