

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)

inputFile="symbol.txt"      #输入文件
geneFile="hub.txt"         #基因列表文件
logFCfilter=1               #logFC临界值
fdrFilter=0.05              #fdr临界值

setwd("D:\\肝细胞癌\\12.hub基因TCGA生存分析\\1.hubGeneExp")     #设置工作目录

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

#正常和肿瘤数目
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
Type=c(rep(1,conNum), rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]


#提取差异的hub基因
row.names(outDiff)=outDiff[,1]
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), as.vector(outDiff[,1]))
immuneDiff=outDiff[sameGene,]
write.table(immuneDiff, file="TCGA.hubDiff.txt", sep="\t", row.names=F, quote=F)

#差异hub基因的表达文件
immuneExp=rbind(ID=colnames(data[sameGene,]), data[sameGene,])
write.table(immuneExp, file="TCGA.hubDiffExp.txt", sep="\t", col.names=F, quote=F)


