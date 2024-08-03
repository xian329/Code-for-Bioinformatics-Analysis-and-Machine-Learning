
rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\ssgsea标志基因集功能分析\\S&R\\3.单基因的相关性检验") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DOSE")
#library(DOSE)

#引用包
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

#输入文件
h.all="ssgseaOut.txt"         #ssgsea富集结果
exp="uniq.symbol.txt"    #肿瘤表达矩阵
hub="hub.txt"     #核心基因文件名称

rt=read.table(exp,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
immune=read.table(h.all, header=T, sep="\t", check.names=F, row.names=1)
immune=t(immune)
data=immune
risk=rt
af=read.table(file = hub,sep = "\t",header = F,check.names = F)
risk=risk[af[,1],]
sameSample=intersect(colnames(risk),rownames(data))
data=data[sameSample,]
risk=risk[,sameSample]
risk=t(risk)
data1=data
outTab=data.frame()
for(immune in colnames(data)){
  for(gene in colnames(risk)){
    x=as.numeric(data[,immune])
    y=as.numeric(risk[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
  }
}
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=8, height=15)
ggplot(outTab, aes(Gene, Immune)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(high="red", mid = "white",low= "green") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +   
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  
        axis.text.y = element_text(size = 8, face = "bold")) +      
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) + 
  scale_x_discrete(position = "bottom")     
dev.off()


