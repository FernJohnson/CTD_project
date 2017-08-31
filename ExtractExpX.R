function(eSet, name, IDs)
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
