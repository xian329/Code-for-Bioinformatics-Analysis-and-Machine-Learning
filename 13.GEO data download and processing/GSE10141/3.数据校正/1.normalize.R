

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#引用包
library(limma)
library(sva)

geoExpFile="GSE10141.txt"      #GEO表达数据文件
setwd("D:\\肝细胞癌\\13.GEO数据下载及处理\\GSE10141\\3.数据校正")      #设置工作目录


#读取geo基因表达文件,并对数据进行处理
rt=read.table(geoExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)

#如果GEO数据没有取log2,会自动对数据取log2
qx=as.numeric(quantile(geo, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo[geo<0]=0
  geo=log2(geo+1)}
geo=normalizeBetweenArrays(geo)

geo[geo<0]=0

#输出矫正后的数据
geoTab=rbind(ID=colnames(geo), geo)
write.table(geoTab,file="GEO.normalize.txt",sep="\t",quote=F,col.names=F)

