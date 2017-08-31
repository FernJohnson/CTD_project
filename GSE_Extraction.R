#Downloading and extracting expression data from GSE GEO entries

library(GEOquery)

#First obtain set of GEO accessions - assuming GDS accessions were found earlier
geogse<-setdiff(georel, geogds)

#There are many more GSEs than GDS; best to deal with them in batches

GSE2<-geogse[11:20]


#Creates esets, and also stores them is a list 'gsestuff'
gsestuff<-c()

#Usually GSEs can be directly downloaded into an eset.
#Or rather, a list containing an eset - needs to extracted, as in b[[k]]
#Some GSEs have multiple esets
#The code below checks for this, and extracts all the esets 

for(i in GSE2[1:10])
{
  b<-getGEO(i, GSEMatrix = TRUE)
  
  if(length(b)>1)
  {
    for(k in length(b))
      gsestuff<-c(gsestuff, b[[k]])
    d<-paste(i,k,"eSet", sep="_")
    assign(d, b[[k]])
  }
  else
  {
    gsestuff<-c(gsestuff, b[[1]])
    c<-paste(i, "eSet", sep="_")
    assign(c, b[[1]])
  }
  
}

#Some accessions won't work.  Can attempt to build eset;

GSE86425<-getGEO(GSE2[7], GSEMatrix = FALSE)


#Procedure outlined fully in the 'Using the GEOquery package' user guide

#First, check all GSMs (samples) have the same platform
gsmlist = Filter(function(gsm) {Meta(gsm)$platform=='GPL11180'},GSMList(GSE86425))
length(gsmlist)

Table(gsmlist[[1]])[1:5,]

Columns(gsmlist[[1]])[1:5,]

probesets <- Table(GPLList(GSE86425)[[1]])$ID

Meta(GSE86425)
phen1<-GSMList(GSE86425)[[1]]
names(phen1@header)

gsmplatforms <- lapply(GSMList(GSE86425),function(x) {Meta(x)$platform})
tail(gsmplatforms)

phenlist<-lapply(gsmlist, function(x) {x@header})
cols<-names(phenlist)
rows<-names(phen1@header)


phenlist$GSM2302279$
  
n<-length(phenlist[[1]])

DF <- structure(phenlist, row.names = c(NA, -n), class = "data.frame")
head(DF)


dim(DF)
rownames(DF)<-rows
colnames(DF)<-cols

data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
{tab <- Table(x)
mymatch <- match(probesets,tab$ID_REF)
return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})

data.matrix[1:5,]

require(Biobase)

rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
exprs(eset2)

annotation(eset2)<-"GPL11180"

GSE86425_eSet<-eset2

#Extracting and processing data
#It is useful to look at the platform information first
#If multiple esets have similar platform information, i.e they have gene symbols,
#they can be batch processed, using a loop at the ExtractExp function

platf<-getGEO(annotation(GSE86425_eSet), AnnotGPL = TRUE)
IDs<-attr(dataTable(platf), "table")[, c("ID", "Gene symbol")]

sets<-c(GSE61728_eSet, GSE34896_eSet, GSE58173_eSet, GSE48280_eSet, GSE46239_eSet,
        GSE39454_eSet, GSE26852_eSet, GSE14997_eSet, GSE11083_eSet)

names<-c("GSE61728", "GSE34896", "GSE58173", "GSE48280", "GSE46239", "GSE39454",
         "GSE26852", "GSE14997", "GSE11083")

for (i in 1:length(sets))
{
  platf<-getGEO(annotation(GSE10144_1_eSet), AnnotGPL = TRUE)
  IDs<-attr(dataTable(platf), "table")[, c("ID", "Symbol")]
  
  f<-ExtractExpX(sets[[i]], names[i], IDs)
  (assign(paste(names[i], "exp", sep = "_"), f))
}


#Then, extract the values, rank transformed values, and phenotypic info
#It is worth doing this manually each time, and checking the table before saving it
#it is much easier to spot problems at this stage

head(GSE61728_exp)

write.csv(GSE61728_exp, file="GSE61728_exp.csv")

GSE61728_exprank<-apply(GSE61728_exp, 2, rank)

head(GSE61728_exprank)

write.csv(GSE61728_exprank, file="GSE61728_exprank.csv")

phen<-t(pData(phenoData(GSE61728_eSet)))

write.csv(phen, file="GSE61728_pheno.csv")

