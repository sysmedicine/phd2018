#### Installs required packages ####
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2)
}
if (!require("PIANO")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("piano")
  library(piano)
}
if (!require(pheatmap)) {   #Corrected
  source("http://bioconductor.org/biocLite.R")
  biocLite("pheatmap")
  library(pheatmap)
}
if (!require("factoextra")) {
  install.packages("factoextra", dependencies = TRUE)
  library(factoextra)
}
if(!require("corrplot")){
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
}
if(!require("Hmisc")){
  install.packages("Hmisc", dependencies = TRUE)
  library(Hmisc)
}
if(!require("ConsensusClusterPlus")){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ConsensusClusterPlus")
}
if(!require("reshape2")){
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}

#### File path ####
#Always run these lines. They are used to identify your file path. Used throughout
if(dir.exists(paste0(getwd(),"/Code/"))){ 
  path=paste0(getwd(),"/Code/")
}else{
  path=paste0(getwd(),"/")
}
pathRaw=paste0(path,'/Raw/')

#### DESeq2 ####
library(DESeq2) #Loads the DESeq2 package to be used.

#Imports counts data
Counts=as.matrix(read.csv(file = paste0(pathRaw,"Liver_HCC_all_counts.txt"), 
                          sep = "\t", row.names = 1))

#Imports metadata
Metadata=read.csv(file = paste0(pathRaw,"Liver_HCC_all_counts_metadata.txt"), 
                  sep = "\t",row.names = 1)

#Shows a preview of what the Counts and Metadata looks like
head(Counts) 
head(Metadata)  #notice that one column is called "type", used below

#Creates an object that includes count and metadata to run DESeq2
deSeqData = DESeqDataSetFromMatrix(countData=Counts, colData=Metadata, design= ~type)

deSeqAnalysis= DESeq(deSeqData) #Performs differential expression. This is the main function

#Extracts information from DESeq analysis considering the comparisons "Cancer" vs "Matched"
res = results(deSeqAnalysis, 
              contrast=c("type","Cancer","Matched"), #labels from contrasts between samples. Specified in Metadata.
              lfcThreshold = 0,altHypothesis="greaterAbs",#Null hypothesis: LFC = 0; Alternative hypothesis: two-sided comparison 
              alpha = 0.01) #FDR considered as threshold for significance

summary(res) #returns a summary of the results

#Outputs DESeq result to a file called "DESeq output.txt". You use this file in PIANO
write.table(as.data.frame(res), file=paste0(path,"DESeq output.txt"),
            sep = "\t",row.names = T,col.names = NA)

#### GSEA - PIANO ####
library(piano)

#Imports DESeq output
DESeqout=read.delim(file=paste0(path,"DESeq output.txt"), row.names = 1, stringsAsFactors = F)
head(DESeqout) #Preview of DESeq output. Notice that all transcripts are annottated with Ensembl ids

#Imports a file used to map Ensembl ids to Gene Symbols
ens2gene=read.delim(paste0(path,"Ensembl2gene.tsv"), row.names = 2, stringsAsFactors = F)

DESeqout[,"Gene"]=ens2gene[row.names(DESeqout),"Gene"] #Adds a column to DESeqout with the gene names

#The following lines prepare the DESeq output to be used in GSEA
#And for simplicity, we take the 1st transcript for genes with multiple Ensembl ids, and drop all duplicates.
DESeqout=DESeqout[!is.na(DESeqout$Gene),] #Excludes Ensembl ids without corresponding gene symbols
DESeqout=DESeqout[!duplicated(DESeqout$Gene),] #A gene symbol may be associated with different Ensembl ids, and some Ensembl ids have no associated gene symbol. 
row.names(DESeqout)=DESeqout$Gene #Assigns the gene symbols as row names of DESeqout 

#Piano only needs LFC and pvalue from DESeq
DESeqout=DESeqout[ ,c('log2FoldChange','pvalue')] 
pval= as.matrix(DESeqout[ ,2]) #extract P as a matrix
fc= as.matrix(DESeqout[ ,1])  #extract fold changes as a matrix
row.names(pval)=row.names(DESeqout)
row.names(fc)=row.names(DESeqout)

# P values and FC identified as NA by DESeq2 are here considered as non-significant
pval[is.na (pval)] <- 1
fc[is.na (fc)] <- 0

#Gene set downloaded from MSigDB: 
# http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C5
# GO biological processes symbols
gset=loadGSC(paste0(path,"c5.bp.v6.2.symbols.gmt")) 

#Main function in PIANO. Uses the pvalues, fold-changes, and gene sets. 
#Assigns FDR as multiple hypothesis correction method
gsaRes <- runGSA(pval,fc,gsc=gset, nPerm = 1000, adjMethod = "fdr")

#Writes the output from PIANO to a file "piano.txt". Setting "save=FALSE" simply shows the result from PIANO
GSAsummaryTable(gsaRes, save=TRUE, file=paste0(path,"piano.txt"))

#### Hierarchical clustering & Heatmap visualization ########
library(pheatmap) #Package to generate heatmaps.

#Reading FPKM data. Column names show conditions (Normal, Cancer)
FPKMdata=as.matrix(read.csv(paste0(pathRaw,"Liver_HCC_all_fpkm_selected.txt"), sep = "\t", row.names = 1))

#Since some of these analyses take time to run, for this test case we will use the 15 genes showing highest variance
varianceData=as.data.frame(apply(FPKMdata,1,var)); #computes variance
colnames(varianceData)="variance" #renames column
topgenes=row.names(varianceData[order(varianceData$variance, decreasing = T),,drop=F])[1:15] #Takes top varying genes
FPKMdata=FPKMdata[topgenes,] #filters table based on those genes

#Computes distances using "euclidean" distance, and the hierarchical clustering using "average" linkage between clusters
#May be used to output the metrics of the hierarchical clustering
hclust(dist(FPKMdata, method = "euclidean"), method = "average") 

# Heatmap of expression data
pheatmap(FPKMdata, 
  scale = "row", #note that data is scaled rowise to show differences between groups
  clustering_distance_rows = "euclidean", #distance method
  clustering_distance_cols = "correlation", #distance method
  clustering_method = "complete" #linkage method
  )

#### PCA ####
library(gplots)
library(RColorBrewer)
library(factoextra)

FPKMdata=as.matrix(read.csv(paste0(pathRaw,"Liver_HCC_all_fpkm_selected.txt"), sep = "\t", row.names = 1))
FPKMdata=FPKMdata[rowMeans(FPKMdata)>1,] #Ignores genes that are undetected in most samples since these do not contribute for data variance

FPKMdata=t(FPKMdata) #Necessary to transpose matrix because columns become features (variables) for PCA

#Adds a column with labels for the replicates
FPKMdata=cbind(FPKMdata,data.frame("group"=c("Normal","Normal","Normal","Normal","Normal","Normal",
                                   "Cancer","Cancer","Cancer","Cancer","Cancer","Cancer","Cancer"))
)
numcols=colnames(FPKMdata)[colnames(FPKMdata)!="group"] #Selects columns with FPKM data

#Main function. Computes the principal components of the numeric data (numcols) 
#Note that data is scale (scale=T) which is very important for PCA
pca=prcomp(FPKMdata[,numcols], scale = T)

#Visualization of PCA
fviz_pca_ind(pca, 
             label="none",#removes sample ids from PCA plot
             addEllipses=TRUE, ellipse.level=0.95, #shows ellipses for 0.95 confidence interval
             habillage=FPKMdata$group,
             title="PCA of all expressed genes"
             )

#### Correlation plots ####
library(corrplot)
library(Hmisc)

#The following lines load the data and process it. We use only 500 genes just to speed up this test case
#Reading FPKM data. Column names show conditions (Normal, Cancer)
FPKMdata=as.matrix(read.csv(paste0(pathRaw,"Liver_HCC_all_fpkm_selected.txt"), sep = "\t", row.names = 1))

#Since some of these analyses take time to run, for this test case we will use the 15 genes showing highest variance
varianceData=as.data.frame(apply(FPKMdata,1,var)); #computes variance
colnames(varianceData)="variance" #renames column
topgenes=row.names(varianceData[order(varianceData$variance, decreasing = T),,drop=F])[1:500] #Takes top varying genes
FPKMdata=FPKMdata[topgenes,] #filters table based on those genes

## These are the main functions, used to build the correlation matrix, 
#with associated correlation coefficient and P
cor.matrix=rcorr(as.matrix(FPKMdata), type="spearman")
cor.matrix.R=cor.matrix$r
cor.matrix.P=cor.matrix$P

#Plots the correlation matrix
corrplot(cor.matrix.R, p.mat = cor.matrix.P, 
         order = "hclust", #columns and rows sorted by hierarchical clustering
         insig = "blank", sig.level = 0.001) #shows correlations with P>0.001 as blanks

#### Consensus clustering ####
library(ConsensusClusterPlus)

Dat=as.matrix(read.csv(paste0(pathRaw,"Liver_HCC_all_fpkm_selected.txt"), sep = "\t", row.names = 1))
res=ConsensusClusterPlus(Dat,maxK=6,pItem=0.8,pFeature=1,
                         clusterAlg="hc",distance="pearson",seed=42,plot="png")


