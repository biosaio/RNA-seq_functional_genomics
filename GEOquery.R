#browseVignettes("SummarizedExperiment")
#browseVignettes("GEOquery")

#GSE122305 
GEOid="GSE122305"          # imposto il mio dataset
gse = GEOquery::getGEO(GEOid,destdir=getwd())[[1]]   # scarico i dati
?getGEO

utils::browseURL(getwd()) #opens the download directory which is based on your project directory

library("SummarizedExperiment")

se = as(gse, "SummarizedExperiment")
se
assays(se)$exprs%>%head
colData(se) 
metadata(se)

getwd() #setwd() to change it
if(GEOid %in% list.files()){"you already have it"}else{getGEOSuppFiles(GEOid)} #get a tar file full of subfiles, look for it in your working directory: 

utils::browseURL(paste(getwd(),GEOid,sep="\\")) #opens the download directory 

#this should work if there is only file and the separator is \t, another common option is a comma
files<-list.files(GEOid,full.names=T)
length(files)
df<-read.table(files[1], skip=0, sep = "\t", header=TRUE, row.names = 1)
df%>%rownames
df%>%head
#this file is same the ftp download from geo omnibus