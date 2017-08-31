#Downloading and extracting expression values from ArrayExpress data 

library(ArrayExpress)

library(oligo)

length(aerel)

#These took quite a bit of trial and error - worth doing them in batches, dependent
#on the R container they are downloaded into

#For example, ExpressionFeatureSets

aestuff<-c()

EFS1<-c("E-MEXP-2681", "E-MEXP-1214", "E-MEXP-2769", "E-MEXP-382")

#First, just download from AE into ExpressionFeatureSet container

for (i in EFS1)
{
  b<-ArrayExpress(i)
  aestuff<-c(aestuff,b)
  a<-gsub("-", "", i)
  assign(paste(a,"obj",sep="_"), b)
}

#It is really worth saving the objects individually - RStudio will be very liable to crashing
#And wiping the environment when you're working with these objects
#This will save downloading them all again

saveRDS(EMEXP382_obj, "EMEXP382_obj.RDS")

EFS1strings<-c("EMEXP2681_eSet", "EMEXP1214_eSet", "EMEXP2769_eSet", "EMEXP382_eSet")

EFSesets<-c()

EMEXP2681_obj<-readRDS("EMEXP2681_obj.RDS")
EMEXP1214_obj<-readRDS("EMEXP1214_obj.RDS")
EMEXP2769_obj<-readRDS("EMEXP2769_obj.RDS")
EMEXP382_obj<-readRDS("EMEXP382_obj.RDS")

aestuff<-c(EMEXP2681_obj, EMEXP1214_obj, EMEXP2769_obj, EMEXP382_obj)

#Use the RMA algorithm to convert the expressionfeaturesets to esets

for (i in 1:length(aestuff))
{
  inorm<-rma(aestuff[[i]])
  EFSesets<-c(EFSesets, inorm)
  assign(EFS1strings[i], inorm)
  
}

#Then, process info 

EFS1plats<-c("GPL96", "GPL570", "GPL570", "GPL96")

#Seeing as there are way more GEO datasets, just look up the platform on there

for(i in 1:length(EFSesets))
{
  

  platf<-getGEO(EFS1plats[i], AnnotGPL = TRUE)
 
  IDs<-attr(dataTable(platf), "table")[, c("ID", "Gene symbol")]
  expr<-exprs(EFSesets[[i]]) 
  
  rows<-rownames(expr)
  
  symb<-rep(0, length(rows))
  for(k in 1:length(rows))
  {
    symb[k]<-as.character(IDs[which(IDs[,1] == rows[k]), 2])
  }
  
  rownames(expr)<-symb
  
  sel<-apply(as.matrix(symb, ncol=1), 1, function(x){if(x=="") return(FALSE) else return(TRUE)})
  
  selExpr<-c()
  selSymb<-c()
  
  for(j in 1:nrow(expr))
  {
    if(sel[j])
    {
      selExpr<-rbind(selExpr, expr[j,])
      selSymb<-rbind(selSymb, symb[j])
    }
  }
  
  rownames(selExpr)<-selSymb
  
  selExprm<-aggregate(selExpr, by=list(row.names(selExpr)), FUN=mean)
  
  genenames<-unique(rownames(selExpr))
  
  genenames<-sort(genenames)
  
  rownames(selExprm)<-genenames
  
  selExprm<-selExprm[,-1]
  
  a<-gsub("-", "", EFS1[i])
  
  assign(paste(a, "exp", sep = "_"), selExprm)
  
}

#Write files, rank data, get pheno data 
write.csv(EMEXP2681_exp, file="EMEXP2681_exp.csv")

EMEXP2681_exprank<-apply(EMEXP2681_exp, 2, rank)

write.csv(EMEXP2681_exprank, file="EMEXP2681_exprank.csv")

phen<-t(pData(phenoData(EMEXP2681_eSet)))

write.csv(phen, file="EMEXP2681_pheno.csv")


#Next, N-Channel sets 

Nchan<-c("E-MTAB-5542", "E-MEXP-3881")

for (i in Nchan)
{
  b<-ArrayExpress(i)
  aestuff<-c(aestuff,b)
  a<-gsub("-", "", i)
  assign(paste(a,"obj",sep="_"), b)
}

saveRDS(EMEXP3881_obj, "EMEXP3881_obj.RDS")
saveRDS(EMTAB5542_obj, "EMTAB5542_obj.RDS")

channelNames(EMTAB5542_obj)

EMTAB5542_E_eSet<-channel(EMTAB5542_obj, "E")
EMTAB5542_Eb_eSet<-channel(EMTAB5542_obj, "Eb")

#Then, can attempt to extract expression data and process
#First eset is downloadable normally 

platf<-getGEO(annotation(EMTAB5542_E_eSet), AnnotGPL = TRUE)
IDs<-attr(dataTable(platf), "table")[, c("ID", "Gene symbol")]

f<-ExtractExpX(EMTAB5542_E_eSet, "EMTAB5542_E", IDs)
assign(paste(names[i], "exp", sep = "_"), f)

#And the files can be transformed, and saved to csv as usual
#The second eset has a problem with the platform
#The platform table was located online and downloaded

cols<-c("ID", "COL", "ROW",	"NAME","SPOT_ID","CONTROL_TYPE","REFSEQ","GB_ACC","LOCUSLINK_ID","GENE_SYMBOL","GENE_NAME","UNIGENE_ID","ENSEMBL_ID","ACCESSION_STRING","CHROMOSOMAL_LOCATION",	"CYTOBAND",	"DESCRIPTION",	"GO_ID", "SEQUENCE")

anno<-read.table("GPL16699-15607.txt", sep="\t", col.names=cols,fill=TRUE)

colnames(anno)<-cols

annot<-anno[-1,]

IDs<-annot[, c("ID", "GENE_SYMBOL")]

expr<-exprs(EMTAB5542_Eb_eSet)


rows<-rownames(expr)

symb<-rep(0, length(rows))
for(k in 1:length(rows))
{
  e<-as.character(IDs[which(IDs[,1] == rows[k]), 2])
  
  if(length(e)>0)
  {symb[k]<-e}
  else
  {symb[k]<-""}
  
}

rownames(expr)<-symb

sel<-apply(as.matrix(symb, ncol=1), 1, function(x){if(x=="" || x==0) return(FALSE) else return(TRUE)})

selExpr<-c()
selSymb<-c()

for(j in 1:nrow(expr))
{
  if(sel[j])
  {
    selExpr<-rbind(selExpr, expr[j,])
    selSymb<-rbind(selSymb, symb[j])
  }
}

rownames(selExpr)<-selSymb

selExprm<-aggregate(selExpr, by=list(row.names(selExpr)), FUN=mean)

genenames<-unique(rownames(selExpr))

genenames<-sort(genenames)

rownames(selExprm)<-genenames

selExprm<-selExprm[,-1]

a<-gsub("-", "", "E-MTAB-5542_Eb")

return(assign(paste(a, "exp", sep = "_"), selExprm))


head(EMTAB5542_Eb_exp)

write.csv(EMTAB5542_Eb_exp, file="EMTAB5542_Eb_exp.csv")

EMTAB5542_Eb_exprank<-apply(EMTAB5542_Eb_exp, 2, rank)

write.csv(EMTAB5542_Eb_exprank, file="EMTAB5542_Eb_exprank.csv")

phen<-t(pData(phenoData(EMTAB5542_Eb_eSet)))

write.csv(phen, file="EMTAB5542_Eb_pheno.csv")

#Additionally, info on attempting to build R objects out of raw or processed
#MAGE-TAB files, if these methods do not work, is given in the guide 'Building
#'Building R objects from ArrayExpress datasets' by Audrey Kaufman.  The code for this
#is not included here, as during this analysis these procedures failed due to missing
#data in the datasets of interest.


