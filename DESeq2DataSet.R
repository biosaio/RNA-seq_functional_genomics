
#example when we have transcripts as rows
gse1<-getGEO("GSE122305",GSEMatrix =TRUE)[[1]]
se1 = as(gse1, "SummarizedExperiment")
metadata(se1)
se1@metadata
se1@colData

getwd()

if("GSE122305" %in% list.files()){"you already have it"}else{getGEOSuppFiles("GSE122305")}
list.files("GSE122305")
utils::browseURL(paste(getwd(),"GSE122305",sep="\\"))

df1<-read.table("GSE122305/GSE122305_Raw_counts.txt.gz", skip=0, sep = "\t", header=TRUE)#, row.names = 1
head(df1)
#your dataframe should contain genes as rows and samples as columns,
#hence we should not use Transcript.ID as rownames, use Gene.Name instead
#but some Gene.Names are duplicated.

#in our case there are no duplicates, hence i comment the section 
#df1$Gene.Name[duplicated(df1$Gene.Name)]
#length(unique(df1$Gene.Name))

#In case we are given trascripts we aggregate by gene_id 

df1%>%head
df1b<-df1[,-1]
rownames(df1b)<-df1[,1]
df1b%>%head
#make df1b%>%colnames and se1@colData%>%rownames match
df1b%>%colnames
se1@colData%>%rownames
rownames(se1@colData)<-colnames(df1b)

# TO ASK PROFESSOR CONFIRMATION ABOUT THIS: 
se1@colData
se1@colData$source_name_ch1 #these are our groups
se1@colData$source_name_ch1<-factor(str_replace_all(se1@colData$source_name_ch1, "[^A-Za-z0-9_]", "_"))
se1@colData$source_name_ch1<-factor(str_replace_all(se1@colData$source_name_ch1, "_skin", ""))

#se1@colData$characteristics_ch1.1<-factor(str_replace_all(se1@colData$characteristics_ch1.1, "[^A-Za-z0-9_]", "."))

#if your dataframe was not integer we round
#in our case no change
#df1b %>% mutate_if(is.numeric,round)%>%->df1b


dds1afresh <- DESeqDataSetFromMatrix(countData = df1b,
                              colData = se1@colData,
                              design = ~ source_name_ch1,
                              metadata=metadata(se1))

#accesing values:
#dds1afresh@metadata
#counts(dds1afresh)
#dds1afresh$source_name_ch1

#dds1bfresh is better design, we will see there is signifant effect due to se1@colData$characteristics_ch1.1
# dds1bfresh <- DESeqDataSetFromMatrix(countData = df1b,
#                                     colData = se1@colData,
#                                     design = ~ characteristics_ch1.1+source_name_ch1,
#                                     metadata=metadata(se1))

#se1@colData$genotype.ch1
#se1@colData$characteristics_ch1.1

#----------------------------
#what if multiple files are downloaded?
# GEOid<-"GSE122305" #example:GSE124252, GSE167118,GSE165361 project1: GSE145445, project2: GSE139419, free: GSE154104<-split into multiple files, GSE159522
# if(GEOid %in% list.files()){"you already have it"}else{getGEOSuppFiles(GEOid)}
# list.files(GEOid,full.names = TRUE)
# files<-list.files(GEOid,full.names = TRUE); files
# untar(files[[1]],list = T)
# ?untar
# untar(files[[1]],exdir=GEOid)
# files<-list.files(GEOid,full.names = TRUE); files
# gse2<-getGEO(GEOid,GSEMatrix =TRUE)[[1]]
# se2= as(gse2, "SummarizedExperiment")
# metadata(se2)[[1]]
# se2@colData
# 
# # \\. just means "." but we need the backslash to say its not a regular expression
# #read file extension
# files
# issplit<-lapply(files, str_split, pattern="\\.", n=2)
# lapply(files,function(x){str_split(x,pattern="\\.", n=2)})
# issplit
# unique(sapply(issplit,function(x){x[[1]][[2]]}))
# 
# list.files(GEOid)
# tsvfiles <- list.files(path=GEOid, pattern = '.txt\\.gz$', ignore.case = TRUE, full.names = TRUE)
# tsvfiles
# 
# 
# tsvlist <- lapply(tsvfiles, read.table, sep = "\t", header=TRUE)
# names(tsvlist)<-tsvfiles   #in this way we can access dataframes by name with tsvlist$...
# #tsv is a list that contains dataframes
# 
# tsvlist%>%lapply(head,n=3)
# 
# #join into a single dataframe
# tsvlist%>%lapply(function(x){x[,"count"]})->listofcols
# listofcols%>%map(head,n=3)   #simile a sapply
# 
# df2<-as.data.frame(listofcols)
# df2%>%head
# colnames(df2)
# 
# se2@colData
# colData(se2)[[1]]
# colData(se2)[[2]]
# colData(se2)[[8]]
# colnames(df2)<-colData(se2)[[1]]
# rownames(df2)<-tsvlist[[1]][,1]
# df2%>%head()
# tsvlist[[1]][,1]
# 
# colData(se2)
# rownames(se2@colData)
# rownames(se2@colData)<-colnames(df2)
# se2@colData$source_name_ch1 #these are our groups
# se2@colData$source_name_ch1<-factor(str_replace_all(se2@colData$source_name_ch1, "[^A-Za-z0-9_]", "."))


#if your dataframe was not integer we round
#in our case no change
#df2 %>% mutate_if(is.numeric,round)->df2


# dds2fresh <- DESeqDataSetFromMatrix(countData = df2,
#                                     colData = se2@colData,
#                                     design = ~ source_name_ch1,
#                                     metadata=metadata(se2))


#dds1afresh, dds1bfresh, dds2fresh
# Save an object to a file

saveRDS(list("source_name_ch1"=dds1afresh), file = "dds.rds")
# Restore the object
ddslist<-readRDS(file = "dds.rds")

ddslist$source_name_ch1@metadata
ddslist$source_name_ch1@colData
counts(ddslist$source_name_ch1)

