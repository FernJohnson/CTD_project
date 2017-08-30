#RNA IDs
#Getting accessions from RNAseq experiments, GEO only
#Code from GEO section of 'IDs.R' used.

library(rentrez)
library(GEOquery)

#Sjogren's syndrome
msearch<-entrez_search(db="gds", term="(Sjogren's Syndrome[MESH]) AND Expression profiling by high throughput sequencing[DataSet Type]", retmax=1000)
pssids<-msearch$ids
length(pssids)

linkedmids<-entrez_link(dbfrom="gds",id=pssids,db='pubmed')

psspmedlinks<-linkedmids$links$gds_pubmed


linkz<-entrez_link(dbfrom='pubmed', id=psspmedlinks, db='gds')

pssuid<-linkz$links$pubmed_gds

#psspmdsum<-entrez_summary(db="gds", id=pssuid)
#psspmed<-extract_from_esummary(psspmdsum, "accession", simplify = TRUE)

psspmed<-psspmedlinks

pssuidsum<-entrez_summary(db="gds", id=pssids)
pssaccsseq<-extract_from_esummary(pssuidsum, "accession", simplify = TRUE)

#Antiphospholipid Syndrome

aplssearch<-entrez_search(db="gds", term="(Antiphospholipid Syndrome[MESH])AND Expression profiling by high throughput sequencing[DataSet Type]", retmax=1000)
aplsids<-aplssearch$ids
length(aplsids)

linkedmids<-entrez_link(dbfrom="gds",id=aplsids,db='pubmed')

aplspmedlinks<-linkedmids$links$gds_pubmed
aplspmed<-aplspmedlinks

linkz<-entrez_link(dbfrom='pubmed', id=aplspmedlinks, db='gds')

aplsuidseq<-linkz$links$pubmed_gds

#aplspmdsum<-entrez_summary(db="gds", id=aplsids)
#aplspmed<-extract_from_esummary(aplspmdsum, "accession", simplify = TRUE)

aplsuidsum<-entrez_summary(db="gds", id=aplsids)

aplsaccsseq<-extract_from_esummary(aplsuidsum, "accession", simplify = TRUE)

#Mixed connective tissue disorders

mctdsearch<-entrez_search(db="gds", term="(Mixed connective tissue disease[MESH])", retmax=1000)
mctdids<-mctdsearch$ids
length(mctdids)

linkedmids<-entrez_link(dbfrom="gds",id=mctdids,db='pubmed')

mctdpmedlinks<-linkedmids$links$gds_pubmed


linkz<-entrez_link(dbfrom='pubmed', id=mctdpmedlinks, db='gds')

mctduid<-linkz$links$pubmed_gds

#mctdpmdsum<-entrez_summary(db="gds", id=mctduid)
#mctdpmed<-extract_from_esummary(mctdpmdsum, "accession", simplify = TRUE)
mctdpmed<-mctdpmedlinks

mctduidsum<-entrez_summary(db="gds", id=mctdids)

mctdsaccsseq<-extract_from_esummary(mctduidsum, "accession", simplify = TRUE)


#Myositis

myosearch<-entrez_search(db="gds", term="(Myositis[MESH]) AND Expression profiling by high throughput sequencing[DataSet Type]", retmax=500)
myoids<-myosearch$ids
length(myoids)

linkedmids<-entrez_link(dbfrom="gds",id=myoids,db='pubmed')

myopmedlinks<-linkedmids$links$gds_pubmed
myopmed<-myopmedlinks

linkz<-entrez_link(dbfrom='pubmed', id=myopmedlinks, db='gds')

myouid<-linkz$links$pubmed_gds

#myopmdsum<-entrez_summary(db="gds", id=myouid)
#myopmed<-extract_from_esummary(myopmdsum, "accession", simplify = TRUE)


myouidsum<-entrez_summary(db="gds", id=myoids)

myoaccsseq<-extract_from_esummary(myouidsum, "accession", simplify = TRUE)

#Systemic Sclerosis

ssysearch<-entrez_search(db="gds", term="(Scleroderma, Systemic[MESH]) AND Expression profiling by high throughput sequencing[DataSet Type] NOT cancer NOT development NOT neural")
ssyids<-ssysearch$ids

linkedmids<-entrez_link(dbfrom="gds",id=ssyids,db='pubmed')

ssypmedlinks<-linkedmids$links$gds_pubmed


linkz<-entrez_link(dbfrom='pubmed', id=ssypmedlinks, db='gds')

ssyuid<-linkz$links$pubmed_gds

#ssypmdsum<-entrez_summary(db="gds", id=ssyuid)
#ssypmed<-extract_from_esummary(ssypmdsum, "accession", simplify = TRUE)
ssypmed<-ssypmedlinks

ssyuidsum<-entrez_summary(db="gds", id=ssyids)

ssyaccsseq<-extract_from_esummary(ssyuidsum, "accession", simplify = TRUE)

#Localised scleroderma

slosearch<-entrez_search(db="gds", term="(Scleroderma, Localized[MESH] AND Expression profiling by high throughput sequencing[DataSet Type]", retmax=1000)
sloids<-slosearch$ids
length(sloids)

linkedmids<-entrez_link(dbfrom="gds",id=sloids,db='pubmed')

slopmedlinks<-linkedmids$links$gds_pubmed


linkz<-entrez_link(dbfrom='pubmed', id=slopmedlinks, db='gds')

slouid<-linkz$links$pubmed_gds

#slopmdsum<-entrez_summary(db="gds", id=slouid)
#slopmed<-extract_from_esummary(slopmdsum, "accession", simplify = TRUE)
slopmed<-slopmedlinks

slouidsum<-entrez_summary(db="gds", id=sloids)

sloaccsseq<-extract_from_esummary(slouidsum, "accession", simplify = TRUE)

#Cutaneous Lupus

slcsearch<-entrez_search(db="gds", term="(Lupus Erythematosus, Cutaneous[MESH] AND Expression profiling by high throughput sequencing[DataSet Type]", retmax=1000)
slcids<-slcsearch$ids

linkedmids<-entrez_link(dbfrom="gds",id=slcids,db='pubmed')

slcpmedlinks<-linkedmids$links$gds_pubmed


linkz<-entrez_link(dbfrom='pubmed', id=slcpmedlinks, db='gds')

slcuid<-linkz$links$pubmed_gds

#slcpmdsum<-entrez_summary(db="gds", id=slcuid)
#slcpmed<-extract_from_esummary(slcpmdsum, "accession", simplify = TRUE)
slcpmed<-slcpmedlinks

slcuidsum<-entrez_summary(db="gds", id=slcids)

slcaccsseq<-extract_from_esummary(slcuidsum, "accession", simplify = TRUE)

#SLE

slesearch<-entrez_search(db="gds", term="(Lupus Erythematosus, Systemic[MESH] AND Expression profiling by high throughput sequencing[DataSet Type] NOT Canis lupus familiaris[organism])", retmax=1000)
sleids<-slesearch$ids

linkedmids<-entrez_link(dbfrom="gds",id=sleids,db='pubmed')

slepmedlinks<-linkedmids$links$gds_pubmed


linkz<-entrez_link(dbfrom='pubmed', id=slepmedlinks, db='gds')

sleuid<-linkz$links$pubmed_gds

slepmed<-slepmedlinks

sleuidsum<-entrez_summary(db="gds", id=sleids)

sleaccsseq<-extract_from_esummary(sleuidsum, "accession", simplify = TRUE)

#only these searches yielded any results

geoaccsseq<-c(sloaccsseq,ssyaccsseq,myoaccsseq, mctdsaccsseq, pssaccsseq, sleaccsseq)

geoaccsseq<-unique(geoaccsseq)


length(geoaccsseq)

unique(geoaccsseq, georel)
