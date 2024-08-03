

rm(list=ls())#清空环境以排除干扰
options(stringsAsFactors = F)#防止文件乱码

#单因素cox分析数据准备

library(limma)

setwd("D:\\肝细胞癌\\12.hub基因TCGA生存分析\\2.exptime") 

#1.提取TCGA中肿瘤样本的表达信息
rt=read.table("TCGA.hubDiffExp.txt",sep = "\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=t(data)
data=avereps(data)

#2.读取临床文件
cli=read.table("survival_data.txt",sep="\t",check.names=F,header=T,row.names=1)     

#3.将survival_data和mRNA_data合并
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="hub_survival.txt",sep="\t",row.names=F,quote=F)


