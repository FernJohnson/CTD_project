#Differential Gene Expression Analysis

#First seperate by dataset source - 1 RNAseq exp, 2 beadchip exps
#WBAL shorter is a cvs file containing sample accessions, and state of each sample - healthy or SLE patient
#This is a subset of data from the full phenotypic set

which(colnames(WBAL_shorter)=="GSM1863749")

RNAphen<-WBAL_shorter[2,1:117]

which(colnames(WBAL_shorter)=="GSM730999")

Beadphen1<-WBAL_shorter[2,118:147]
Beadphen2<-WBAL_shorter[2,148:268]

be1<-colnames(Beadphen1)
be2<-colnames(Beadphen2)
rn<-colnames(RNAphen)


#The first beadchip set

WBALB1_exp<-fullexp[,be1]
amatb1<-data.matrix(WBALB1_exp, rownames.force = TRUE)

#Remove any rows with no values

amatb12<-na.omit(amatb1)


#Create design matrix 
designb1<-cbind(healthy=rep(0,30),SLE=rep(0,30))
rownames(designb1)<-be1

for (i in 1:length(be1))
{
  if(Beadphen1[i]=="SLE")
  {
    designb1[i,2]<-1
  }
  else
  {
    designb1[i,1]<-1
  }
}

#Rerank gene expression values
amatb12<-apply(amatb12, 2, rank)

#Calculate ArrayWeights
awb1<-arrayWeights(amatb12, design = designb1)
barplot(awb1, xlab="Array", ylab="weight", col="white")
abline(h=1, lwd=1, lty=2) 


#Fit linear model
fitwb1<-lmFit(amatb12, weights = awb1, design = designb1)
fitwb1<-eBayes(fitwb1)

cont.matrix<-makeContrasts(healthyvsSLE=SLE-healthy, levels=designb1)
fitwb12<-contrasts.fit(fitwb1, cont.matrix)

#get topgenes

topgenesb1<-topTable(fitwb12, number=100,adjust="BH")


f_amatb12<-amatb12[rownames(topgenesb1),]


#####

#Now, one of the RNA datasets
#Same analysis, just different data

WBALR_exp<-fullexp[,rn]
amatr<-data.matrix(WBALR_exp, rownames.force = TRUE)

rowSums(is.na(amatr))

amatr2<-na.omit(amatr)

designr<-cbind(healthy=rep(0,117),SLE=rep(0,117))
rownames(designr)<-rn

for (i in 1:length(rn))
{
  if(RNAphen[i]=="SLE")
  {
    designr[i,2]<-1
  }
  else
  {
    designr[i,1]<-1
  }
}

amatr2<-apply(amatr2, 2, rank)

awr<-arrayWeights(amatr2, design = designr)
barplot(awr, xlab="Array", ylab="weight", col="white")
abline(h=1, lwd=1, lty=2) 

fitwr<-lmFit(amatr2, weights = awr, design = designr)
fitwr<-eBayes(fitwr)

cont.matrix<-makeContrasts(healthyvsSLE=SLE-healthy, levels=designr)
fitwr2<-contrasts.fit(fitwr, cont.matrix)

topgenesr<-topTable(fitwr2, number=100,adjust="BH")

write.table(topgenesb1, file="genelistBEAD1.csv", sep=',',row.names=TRUE, col.names=TRUE, quote=FALSE)


f_amatr2<-amatr2[rownames(topgenesr),]

#Seperate 'healthy' and 'SLE' genes based on log fold change

healthygenes<-c()
SLEgenes<-c()

for(i in rownames(topgenesr))
{
  val<-topgenesr[i,1]
  if(val<0)
  {
    healthygenes<-c(healthygenes, i)
  }
  else
  {
    SLEgenes<-c(SLEgenes, i)
  }
}

#Write lists of gene names, ready for PANTHER analysis

write.table(healthygenes, file="genelistrhealthy.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(SLEgenes, file="genelistrSLE.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


#Analysis on combined datasets

WBAL_exp<-fullexp[,colnames(WBAL_shorter)]
amat1<-data.matrix(WBAL_exp, rownames.force = TRUE)

amat2<-na.omit(amat1)

design<-cbind(healthy=rep(0,268),SLE=rep(0,268))

rownames(design)<-colnames(WBAL_shorter)

for (i in 1:length(a))
{
  if(WBAL_shorter[i]=="SLE")
  {
    design[i,2]<-1
  }
  else
  {
    design[i,1]<-1
  }
}

amat3<-apply(amat2, 2, rank)


aw<-arrayWeights(amat3, design = design)
barplot(aw, xlab="Array", ylab="weight", col="white")
abline(h=1, lwd=1, lty=2) 


fitw2<-lmFit(amat3, weights = aw, design = design)
fitw2<-eBayes(fitw2)

cont.matrix<-makeContrasts(healthyvsSLE=SLE-healthy, levels=design)
fitw3<-contrasts.fit(fitw2, cont.matrix)

topgenes<-topTable(fitw3, number=100, adjust="BH")
rownames(topgenes)

f_amat3<-amat3[rownames(topgenes),]

#Cluster samples based on expression of genes identified by topGenes

colv<-as.dendrogram(hclust(as.dist(1-cor(f_amat3)), method="complete"))


library(dendextend)

colhorder<-labels(colv)

zeros<-rep(0, 268)
liness<-rep("|", 268)

colv %>% labels

colv %>% set("labels",zeros) %>%
  set("labels_col", colcouls) %>%  
  plot(main="colors")

colhorder

first<-colhorder[1:55]
second<-colhorder[56:151]
third<-colhorder[152:268]

colv %>% set("labels", colhorder)

colv %>% set('by_labels_branches_col', value=first, TF_values=c(3,Inf)) %>%
  set('by_labels_branches_col', value=second, TF_values=c(8,Inf)) %>%
  set('by_labels_branches_col', value=third, TF_values=c(6, Inf)) %>%
  rotate(ncolhorder) %>%
  set("labels",liness) %>%
  set("labels_cex", 3) %>%
  set("labels_col", ncolcouls) %>%  
  plot(main="Dendrogram of Gene Expression Samples")

colored_bars(colors=colcouls, dend = colv)
legend("topright", leglabs, cex = 0.6,fill = c(3, 8, 6))
legend(1000,100, c("SLE", "Healthy"), cex=0.6, fill=c(2, 4))

leglabs<-c("GSE29536", "GSE22098", "GSE72509")

