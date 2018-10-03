#### Required packages ####
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
if (!require(pheatmap)) {
  install.packages(pheatmap)
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

#### Other functions ####
p=paste0
path=p(getwd(),'/Code/')
pathRaw=p(path,'Raw/')

#### Heatmaps of fold-change, LFC, t-tests, and adjusting p-values for multiple hypothesis correction ####
library(gplots)
library(RColorBrewer)

din=read.csv(p(path,"Book1.txt"), sep = "\t", row.names = 1)
norm=c('N1','N2','N3','N4')
dis=c('D1','D2','D3','D4')

din[,'N_av']=signif(rowMeans(din[,norm]),2)
din[,'D_av']=signif(rowMeans(din[,dis]),2)

heatmap.2(as.matrix(din),
          scale="row",
          Colv=NA, Rowv = NA,dendrogram = 'none',
          trace="none",
          cellnote = din,
          notecex=1.0,
          notecol="black"
        )

din[,'FC']=signif(din$D_av/din$N_av,2)
din[,'LFC']=signif(log2(din$FC),2)
heatmap.2(as.matrix(din[,c('FC','LFC')]),
          scale="column",
          Colv=NA, Rowv = NA,dendrogram = 'none',
          trace="none",
          cellnote = as.matrix(din[,c('FC','LFC')]),
          notecex=1.0,
          notecol="black"
)

myP=c()
for(rr in row.names(din)){
  myP=c(myP,t.test(din[rr,norm], din[rr,dis])$p.value)
}
din[,"P"]=signif(myP,2)
din[,"FDR"]=signif(p.adjust(din$P),2)
heatmap.2(as.matrix(din[,c('P','FDR')]),
          scale="column",
          Colv=NA, Rowv = NA,dendrogram = 'none',
          trace="none",
          cellnote = as.matrix(din[,c('P','FDR')]),
          notecex=1.0,
          notecol="black"
)

# write.table(din,p(path,"temp.txt"),sep = "\t")

#### DESeq2 ####
library(DESeq2)


#Imports and prepares data
Dat=as.matrix(read.csv(p(pathRaw,"Liver_HCC_all_counts_selected.txt"), sep = "\t", row.names = 1))
metaDat=t(Dat[1,,drop=F])
countDat=Dat[2:nrow(Dat),]
class(countDat)="numeric"

#Creates a DESeq subclass, specifying the data and metadata
deSeqData = DESeqDataSetFromMatrix(countData=countDat, colData=metaDat, design= ~type)
deSeqData= DESeq(deSeqData) #Performs differential expression. This is the main function
res = results(deSeqData,#Extracts information from DESeq analysis
              contrast=c("type","Cancer","Matched"),
              lfcThreshold = 0,altHypothesis="greaterAbs",
              alpha = 0.05)
summary(res)

write.table(as.data.frame(res), file=p(path,"DESeq output.txt"),#Outputs DESeq result
            sep = "\t",row.names = T,col.names = NA)

#### GSEA - PIANO ####
library(piano)

din=read.delim(file=p(path,"DESeq output.txt"), row.names = 1, stringsAsFactors = F)
ens2gene=read.delim(p(path,"Ensembl2gene.tsv"), row.names = 2, stringsAsFactors = F)
din[,"Gene"]=ens2gene[row.names(din),"Gene"]

#Different ensembl ids may pertain to the same gene. 
#In this test case and for simplicity, we simply drop duplicates and take the 1st one.
din=din[!is.na(din$Gene),]
din=din[!duplicated(din$Gene),]
row.names(din)=din$Gene
din=din[ ,c('log2FoldChange','pvalue')]
pval= din[ ,2] #extract p
fc= din[ ,1]  #extract fold changes
pval = as.matrix(pval); row.names(pval)=row.names(din)
fc = as.matrix(fc); row.names(fc)=row.names(din)
pval[is.na (pval)] <- 1
fc[is.na (fc)] <- 0

gset=loadGSC(p(path,"c5.bp.v6.2.symbols.gmt")) 

gsaRes <- runGSA(pval,fc,gsc=gset, nPerm=1000)

GSAsummaryTable(gsaRes, save=TRUE, file=p(path,"piano.txt"))

#### Hierarchical clustering ########
library(pheatmap)
Dat=as.matrix(read.csv(p(pathRaw,"Liver_HCC_all_fpkm_selected.txt"), sep = "\t", row.names = 1))
myvar=as.data.frame(apply(Dat,1,var)); colnames(myvar)="variance"
myvar=row.names(myvar[order(myvar$variance, decreasing = T),,drop=F])[1:15] #Takes top varying genes

#Computes distances using "euclidean" distance, 
#and the hierarchical clustering using "average" linkage between clusters
hclust(dist(Dat, method = "euclidean"), method = "average") 

pheatmap(
  Dat[myvar,], scale = "row",
  clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "correlation",
  clustering_method = "complete"
  )

#### PCA ####
library(gplots)
library(RColorBrewer)
library(factoextra)

Dat=as.matrix(read.csv(p(pathRaw,"Liver_HCC_all_fpkm_selected.txt"), sep = "\t", row.names = 1))
Dat=Dat[rowMeans(Dat)>1,] #Ignores genes that are undetected since these do not contribute for data variance

#Transposes matrix and adds label column
Dat=t(Dat) #Necessary because columns become features (variables) for PCA
Dat=cbind(Dat,data.frame("group"=c("Normal","Normal","Normal","Normal","Normal","Normal",
                                   "Cancer","Cancer","Cancer","Cancer","Cancer","Cancer","Cancer"))
)
numcols=colnames(Dat)[colnames(Dat)!="group"] #Only numeric cols

pca=prcomp(Dat[,numcols], scale = T)
fviz_pca_ind(pca, 
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95,
             habillage=Dat$group,
             title="PCA of all expressed genes"
             )

#### Correlation plots ####
library(corrplot)
library(Hmisc)

Dat=as.matrix(read.csv(p(pathRaw,"Liver_HCC_all_fpkm_selected.txt"), sep = "\t", row.names = 1))
Dat=Dat[rowMeans(Dat)>1,] #Ignores genes that are undetected since these do not contribute for data variance
myvar=as.data.frame(apply(Dat,1,var)); colnames(myvar)="variance"
myvar=row.names(myvar[order(myvar$variance, decreasing = T),,drop=F])[1:500] #Takes top varying genes

cor.matrix=rcorr(as.matrix(Dat[myvar,]), type="spearman")
cor.matrix.R=cor.matrix$r
cor.matrix.P=cor.matrix$P

corrplot(cor.matrix.R, p.mat = cor.matrix.P, order = "hclust",
         insig = "blank", sig.level = 0.001)

