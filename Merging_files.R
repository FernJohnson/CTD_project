#Merging tables of values together

#Put all the tables in a file together - put the path to this in 'mypath'

filenames=list.files(path=mypath, full.names=TRUE)
datalist= lapply(filenames, function(x){read.table(file=x,header=T, row.names=1, sep = ",")})
l<-(length(datalist)-1)

v<-datalist[[1]]

#Sometimes the same samples might appear in multiple tables, so these will be removed
#by the 'setdiff' function

#After each merge a new 'row.names' column is produced, this has to be removed

for (i in 1:l)   #(length(datalist)-1))
{
  
  u<-datalist[[i+1]]
  v<-merge(v, u[setdiff(colnames(u),colnames(v))], by=0, all=TRUE)
  rownames(v)<-v$Row.names
  v$Row.names<-NULL

}