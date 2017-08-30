#Gaining CTD database accessions

#GEO
#Getting accessions from microarray experiments only

library(rentrez)
library(GEOquery)

#Sjogren's syndrome

#Simply returns NCBI IDs from search using appropriate MESH term
msearch<-entrez_search(db="gds", term="(Sjogren's Syndrome[MESH]) AND expression profiling by array[DataSet Type]", retmax=1000)
pssids<-msearch$ids
length(pssids)

#Code for extracting pubmed accessions associated with GEO entries

#linkedmids<-entrez_link(dbfrom="gds",id=pssids,db='pubmed')

#psspmedlinks<-linkedmids$links$gds_pubmed

#Code for extracting GEO accessions for entries ONLY associated with a pubmed ID (not used)

#linkz<-entrez_link(dbfrom='pubmed', id=psspmedlinks, db='gds')

#pssuid<-linkz$links$pubmed_gds

#psspmdsum<-entrez_summary(db="gds", id=pssuid)
#psspmed<-extract_from_esummary(psspmdsum, "accession", simplify = TRUE)

#psspmed<-psspmedlinks

#Extract the GEO accession from the search summary, to access GEO entries later
pssuidsum<-entrez_summary(db="gds", id=pssids)
pssaccs<-extract_from_esummary(pssuidsum, "accession", simplify = TRUE)

#Antiphospholipid Syndrome

aplssearch<-entrez_search(db="gds", term="(Antiphospholipid Syndrome[MESH])AND expression profiling by array[DataSet Type]", retmax=1000)
aplsids<-aplssearch$ids
length(aplsids)

#linkedmids<-entrez_link(dbfrom="gds",id=aplsids,db='pubmed')

#aplspmedlinks<-linkedmids$links$gds_pubmed
#aplspmed<-aplspmedlinks

#linkz<-entrez_link(dbfrom='pubmed', id=aplspmedlinks, db='gds')

#aplsuid<-linkz$links$pubmed_gds

#aplspmdsum<-entrez_summary(db="gds", id=aplsids)
#aplspmed<-extract_from_esummary(aplspmdsum, "accession", simplify = TRUE)

aplsuidsum<-entrez_summary(db="gds", id=aplsids)

aplsaccs<-extract_from_esummary(aplsuidsum, "accession", simplify = TRUE)

#Mixed connective tissue disorders

mctdsearch<-entrez_search(db="gds", term="(Mixed connective tissue disease[MESH])", retmax=1000)
mctdids<-mctdsearch$ids
length(mctdids)

#linkedmids<-entrez_link(dbfrom="gds",id=mctdids,db='pubmed')

#mctdpmedlinks<-linkedmids$links$gds_pubmed


#linkz<-entrez_link(dbfrom='pubmed', id=mctdpmedlinks, db='gds')

#mctduid<-linkz$links$pubmed_gds

#mctdpmdsum<-entrez_summary(db="gds", id=mctduid)
#mctdpmed<-extract_from_esummary(mctdpmdsum, "accession", simplify = TRUE)
#mctdpmed<-mctdpmedlinks

mctduidsum<-entrez_summary(db="gds", id=mctdids)

mctdsaccs<-extract_from_esummary(mctduidsum, "accession", simplify = TRUE)


#Myositis

myosearch<-entrez_search(db="gds", term="(Myositis[MESH]) AND expression profiling by array[DataSet Type]", retmax=500)
myoids<-myosearch$ids
length(myoids)

#linkedmids<-entrez_link(dbfrom="gds",id=myoids,db='pubmed')

#myopmedlinks<-linkedmids$links$gds_pubmed
#myopmed<-myopmedlinks

#linkz<-entrez_link(dbfrom='pubmed', id=myopmedlinks, db='gds')

#myouid<-linkz$links$pubmed_gds

#myopmdsum<-entrez_summary(db="gds", id=myouid)
#myopmed<-extract_from_esummary(myopmdsum, "accession", simplify = TRUE)


myouidsum<-entrez_summary(db="gds", id=myoids)

myoaccs<-extract_from_esummary(myouidsum, "accession", simplify = TRUE)

#Systemic Sclerosis

ssysearch<-entrez_search(db="gds", term="(Scleroderma, Systemic[MESH]) AND expression profiling by array[DataSet Type] NOT cancer NOT development NOT neural")
ssyids<-ssysearch$ids

#linkedmids<-entrez_link(dbfrom="gds",id=ssyids,db='pubmed')

#ssypmedlinks<-linkedmids$links$gds_pubmed


#linkz<-entrez_link(dbfrom='pubmed', id=ssypmedlinks, db='gds')

#ssyuid<-linkz$links$pubmed_gds

#ssypmdsum<-entrez_summary(db="gds", id=ssyuid)
#ssypmed<-extract_from_esummary(ssypmdsum, "accession", simplify = TRUE)
#ssypmed<-ssypmedlinks

ssyuidsum<-entrez_summary(db="gds", id=ssyids)

ssyaccs<-extract_from_esummary(ssyuidsum, "accession", simplify = TRUE)

#Localised scleroderma

slosearch<-entrez_search(db="gds", term="(Scleroderma, Localized[MESH] AND expression profiling by array[DataSet Type]", retmax=1000)
sloids<-slosearch$ids
length(sloids)

#linkedmids<-entrez_link(dbfrom="gds",id=sloids,db='pubmed')

#slopmedlinks<-linkedmids$links$gds_pubmed


#linkz<-entrez_link(dbfrom='pubmed', id=slopmedlinks, db='gds')

#slouid<-linkz$links$pubmed_gds

#slopmdsum<-entrez_summary(db="gds", id=slouid)
#slopmed<-extract_from_esummary(slopmdsum, "accession", simplify = TRUE)
slopmed<-slopmedlinks

slouidsum<-entrez_summary(db="gds", id=sloids)

sloaccs<-extract_from_esummary(slouidsum, "accession", simplify = TRUE)

#Cutaneous Lupus

slcsearch<-entrez_search(db="gds", term="(Lupus Erythematosus, Cutaneous[MESH] AND expression profiling by array[DataSet Type]", retmax=1000)
slcids<-slcsearch$ids

#linkedmids<-entrez_link(dbfrom="gds",id=slcids,db='pubmed')

#slcpmedlinks<-linkedmids$links$gds_pubmed


#linkz<-entrez_link(dbfrom='pubmed', id=slcpmedlinks, db='gds')

#slcuid<-linkz$links$pubmed_gds

#slcpmdsum<-entrez_summary(db="gds", id=slcuid)
#slcpmed<-extract_from_esummary(slcpmdsum, "accession", simplify = TRUE)
slcpmed<-slcpmedlinks

slcuidsum<-entrez_summary(db="gds", id=slcids)

slcaccs<-extract_from_esummary(slcuidsum, "accession", simplify = TRUE)

#SLE

slesearch<-entrez_search(db="gds", term="(Lupus Erythematosus, Systemic[MESH] AND expression profiling by array[DataSet Type] NOT Canis lupus familiaris[organism])", retmax=1000)
sleids<-slesearch$ids

#linkedmids<-entrez_link(dbfrom="gds",id=sleids,db='pubmed')

#slepmedlinks<-linkedmids$links$gds_pubmed


#linkz<-entrez_link(dbfrom='pubmed', id=slepmedlinks, db='gds')

#sleuid<-linkz$links$pubmed_gds

#slepmdsum<-entrez_summary(db="gds", id=sleuid)
#slepmed<-extract_from_esummary(slepmdsum, "accession", simplify = TRUE)

#slepmed<-slepmedlinks

sleuidsum<-entrez_summary(db="gds", id=sleids)

sleaccs<-extract_from_esummary(sleuidsum, "accession", simplify = TRUE)



#Combine all the accessions from each search (and pubmed IDs, if desired)

geopmed<-c(slepmed, slcpmed, slopmed, ssypmed, myopmed, mctdpmed, aplspmed, psspmed)
geoaccs<-c(sleaccs,slcaccs,sloaccs,ssyaccs,myoaccs, mctdsaccs, aplsaccs, pssaccs)


#Extract a list of unique accessions and IDs
geoaccs<-unique(geoaccs)
geopmed<-unique(geopmed)

#Check the length
length(geoaccs)
length(geopmed)



#ArrayExpress 

library(ArrayExpress)

#Sjogren's Syndrome

#Search of EPO term returns ArrayExpress ID and pubmed ID (if applicable)
sssets=queryAE(keywords = "Sjogren syndrome")

sspmid<-sssets$PubmedID
ssID<-sssets$ID


#Code for extracting accessions only with associated pubmed IDs (not used)

#sspmids<-as.character(sspmid)
#ssIDs<-as.character(ssID)
#pm<-c()
#ae<-c()



#for (i in 1:length(ssIDs))
#{
#  if(sspmids[i]!="NA")
#    {
#    ae<-c(ae, ssIDs[i])
#    pm<-c(pm, sspmids[i])
#  }
#}

#ssaepm<-cbind(pm, ae)

#Lupus nephritis

lnsets=queryAE(keywords = "lupus nephritis")

lnpmid<-lnsets$PubmedID
lnID<-lnsets$ID

#lnpmids<-as.character(lnpmid)
#lnIDs<-as.character(lnID)
#pm<-c()
#ae<-c()



# for (i in 1:length(lnIDs))
# {
#   if(lnpmids[i]!="NA")
#   {
#     ae<-c(ae, lnIDs[i])
#     pm<-c(pm, lnpmids[i])
#   }
# }
# 
# lnaepm<-cbind(pm, ae)

#SLE

slesets=queryAE(keywords = "systemic lupus erythematosus")

slepmid<-slesets$PubmedID
sleID<-slesets$ID

# slepmids<-as.character(slepmid)
# sleIDs<-as.character(sleID)
# pm<-c()
# ae<-c()
# 
# 
# 
# for (i in 1:length(sleIDs))
# {
#   if(slepmids[i]!="NA")
#   {
#     ae<-c(ae, sleIDs[i])
#     pm<-c(pm, slepmids[i])
#   }
# }
# 
# sleaepm<-cbind(pm, ae)

#systemic scleroderma



syssets=queryAE(keywords = "systemic scleroderma")

syspmid<-syssets$PubmedID
sysID<-syssets$ID

# syspmids<-as.character(syspmid)
# sysIDs<-as.character(sysID)
# pm<-c()
# ae<-c()
# 
# 
# 
# for (i in 1:length(sysIDs))
# {
#   if(syspmids[i]!="NA")
#   {
#     ae<-c(ae, sysIDs[i])
#     pm<-c(pm, syspmids[i])
#   }
# }
# 
sysaepm<-cbind(pm, ae)

#scleroderma (local)

ssets=queryAE(keywords = "scleroderma")

spmid<-ssets$PubmedID
sID<-ssets$ID

# spmids<-as.character(spmid)
# sIDs<-as.character(sID)
# pm<-c()
# ae<-c()
# 
# 
# 
# for (i in 1:length(sIDs))
# {
#   if(spmids[i]!="NA")
#   {
#     ae<-c(ae, sIDs[i])
#     pm<-c(pm, spmids[i])
#   }
# }
# 
# saepm<-cbind(pm, ae)

#cutaneous lupus

clsets=queryAE(keywords = "cutaneous lupus erythematosus")

clpmid<-clsets$PubmedID
clID<-clsets$ID

# clpmids<-as.character(clpmid)
# clIDs<-as.character(clID)
# pm<-c()
# ae<-c()
# 
# 
# 
# for (i in 1:length(clIDs))
# {
#   if(clpmids[i]!="NA")
#   {
#     ae<-c(ae, clIDs[i])
#     pm<-c(pm, clpmids[i])
#   }
# }
# 
# claepm<-cbind(pm, ae)

#myositis

myosets=queryAE(keywords = "myositis")

myopmid<-myosets$PubmedID
myoID<-myosets$ID

# myopmids<-as.character(myopmid)
# myoIDs<-as.character(myoID)
# pm<-c()
# ae<-c()
# 
# 
# 
# for (i in 1:length(myoIDs))
# {
#   if(myopmids[i]!="NA")
#   {
#     ae<-c(ae, myoIDs[i])
#     pm<-c(pm, myopmids[i])
#   }
# }
# 
# myoaepm<-cbind(pm, ae)

#antiphospholipid syndrome

apssets=queryAE(keywords = "primary antiphosolipid syndrome")

apspmid<-apssets$PubmedID
apsID<-apssets$ID

# apspmids<-as.character(apspmid)
# apsIDs<-as.character(apsID)
# pm<-c()
# ae<-c()
# 
# 
# 
# for (i in 1:length(apsIDs))
# {
#   if(apspmids[i]!="NA")
#   {
#     ae<-c(ae, apsIDs[i])
#     pm<-c(pm, apspmids[i])
#   }
# }
# 
# apsaepm<-cbind(pm, ae)

#mixed connective tissue disease

mixsets=queryAE(keywords = "mixed connective tissue disease")

mixpmid<-mixsets$PubmedID
mixID<-mixsets$ID

# mixpmids<-as.character(mixpmid)
# mixIDs<-as.character(mixID)
# pm<-c()
# ae<-c()
# 
# 
# 
# for (i in 1:length(mixIDs))
# {
#   if(mixpmids[i]!="NA")
#   {
#     ae<-c(ae, mixIDs[i])
#     pm<-c(pm, mixpmids[i])
#   }
# }
# 
# mixaepm<-cbind(pm, ae)

#Combining the lists of AE entries with associated pubmed entries (unused)

#aepubmed<-rbind(mixaepm, apsaepm, myoaepm, claepm, saepm, sysaepm, sleaepm, lnaepm, ssaepm)






#Combine lists of AE accessions, extract unique set

aeIDs<-c(mixIDs, apsIDs, myoIDs, clIDs, sIDs, sysIDs, sleIDs, lnIDs, ssIDs)
aeIDs<-unique(aeIDs)

#Remove entries imported from the GEO

aeonly<-c()
for(i in aeIDs)
{
  
  
  if(grepl("GEO", i, fixed = TRUE)==FALSE)
    aeonly<-c(aeonly,i)
  
  
  
}

#aeonly

#After manually checking id entries, remove innapropriate records from GEO and AE
#These include data with a non-microarray experiment, a non-mammalian organism,
#or not CTD related.

removeae<-c("E-TABM-185", "E-MTAB-3732", "E-MTAB-62", "E-MTAB-4475", "E-MTAB-1788", "E-MTAB-4441", "E-ERAD-179")
aerel<-setdiff(aeonly, removeae)

removegeo<-c("GSE50892", "GSE34154", "GSE58613", "GSE6520", "GDS3785", "GSE68766", "GSE45991",
             "GSE30119", "GSE26145", "GSE16129", "GSE9340", "GSE19151", "GSE77599")

georel<-setdiff(geoaccs, removegeo)
length(georel)
