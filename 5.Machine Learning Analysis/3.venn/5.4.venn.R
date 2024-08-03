rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\疾病特征基因\\4.venn") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DOSE")
#library(DOSE)

library(randomcoloR)
library(venn) 

#文件名
A="LASSO"
B="RandomForest"
C="SVM-REF"

#构建一个列表
geneList=list()
rt=read.table(paste0(A,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[A]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("1",uniqLength,sep=" "))

rt=read.table(paste0(B,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[B]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))

rt=read.table(paste0(C,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[C]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))

mycol <- distinctColorPalette(3)
pdf(file="hub.pdf",width=5,height=5)                                                
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)          
write.table(file="hub.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F) 

