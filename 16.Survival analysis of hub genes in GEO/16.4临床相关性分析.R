##########临床相关性分析
#本次挑选的单基因为SMG9，DDX49，NR2C2AP，BSG，OIT3，TPM2
####以SMG9为例
#数据处理
library(limma)
library(ComplexHeatmap)
library(ggpubr)
library(tidyverse)
exp4<-read.table('symbol.txt',sep = '\t',header = T)
exp4=exp4[!duplicated(exp4$id),]
rownames(exp4)<-exp4$id
gene<-exp4[exp4$id=='SMG9',]
gene=gene[,-1]
gene=as.data.frame(t(gene))
#读取表型数据
phe<-read.delim('./clinical.tsv',header = T,sep='\t' )
phe1=phe[,c(2,26:29)]
names(phe1)<-c('id','M','N','stage','T')
table(phe1$M)
#N型数据合并
table(phe1$N)
phe1$N=ifelse(phe1$N=="'--",NA,phe1$N)
table(phe1$T)
phe1$T=ifelse(phe1$T=="'--",NA,phe1$T)
phe1$T=gsub("[ab]$","",phe1$T)
table(phe1$stage)
phe1$stage=ifelse(phe1$stage=="'--",NA,phe1$stage)
#采用正则表达将字母进行删除
phe1$stage=gsub("[ABC]$","",phe1$stage)
phe1$stage=gsub('Stage ','',phe1$stage)
phe1=na.omit(phe1)
gene$sample<-rownames(gene)
gene$sample<-str_replace_all(gene$sample,"\\.","-")
gene$sample=substring(gene$sample,1,12)
phe2=merge(gene,phe1,by.x='sample',by.y = 'id')
phe2=unique(phe2)
library(dplyr)
phe2<-phe2%>%group_by(sample)%>%slice_head(n=1)%>%ungroup()  #选择第一个样本，丢弃第二个
phe2=as.data.frame(phe2)
rownames(phe2)<-phe2$sample
phe2=phe2[,-1]
phe2$SMG9=as.numeric(phe2$SMG9)
gene=colnames(phe2)[1]
tumorData=as.matrix(phe2[gene])
data=avereps(tumorData)
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
Type=factor(Type, levels=c("Low","High"))
data1=cbind(as.data.frame(data), Type)
data1=data1[order(data1[,gene]),] 

samSample=intersect(row.names(data), row.names(phe2))
data=data[samSample,,drop=F]
phe2=phe2[samSample,,drop=F]
rt1=cbind(data, phe2)
rt1=rt1[,-2]
rt1=na.omit(rt1)
for(clinical in colnames(rt1[,2:ncol(rt1)])){
  data=rt1[c(gene, clinical)]
  colnames(data)=c(gene, "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  
  group1=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group1)
  comp=combn(group1,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  boxplot=ggboxplot(data, x="clinical", y=gene, fill="clinical",
                    xlab=clinical,
                    ylab=paste(gene, " expression"),
                    legend.title=clinical)+ 
    stat_compare_means(comparisons = my_comparisons)
  
  pdf(file=paste0(gene,"_clinicalCor_", clinical, ".pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}
















