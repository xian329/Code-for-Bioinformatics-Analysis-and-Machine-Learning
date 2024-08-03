

rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\7.ssgsea标志基因集功能分析\\1.计算得分") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DOSE")
#library(DOSE)

library(GSVA)
library(limma)
library(GSEABase)
expFile="symbol.txt"       #表达矩阵
gmtFile="immune.gmt"       #标志基因集gmt文件    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=normalizeBetweenArrays(mat)
mat=mat[rowMeans(mat)>0,]

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut, file="ssgseaOut.txt", sep="\t", quote=F, col.names=F)

