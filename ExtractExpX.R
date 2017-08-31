#A function that takes in an eSet, string (for naming the table of values), and the IDs, a table of
#Probesets and corresponding gene symbols.  Returns gene expression values annotated with gene
#symbols (or other identifier specfied), with multiple observations aggregated by their mean. 

ExtractExp <- function(eSet, name, IDs)
{
  
  
  expr<-exprs(eSet)
  
  name1<-name
  
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
  
  return(assign(paste(name1, "exp", sep = "_"), selExprm))
  
}

#the function was compiled to improve it's speed

library(compiler)

ExtractExpX<-cmpfun(ExtractExp)
