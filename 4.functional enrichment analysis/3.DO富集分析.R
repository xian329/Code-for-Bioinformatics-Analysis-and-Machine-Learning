rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\富集分析") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("DOSE")

pvalueFilter=0.05        
qvalueFilter=1       
showNum=20

rt=read.table("disease.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
kk <- enrichDO(gene =gene, ont = "DO",pvalueCutoff = 1, qvalueCutoff = 1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="DO.xls",sep="\t",quote=F,row.names = F)
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}
pdf(file="DO_barplot.pdf",width =9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()
pdf(file="DO_bubble.pdf",width =9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

