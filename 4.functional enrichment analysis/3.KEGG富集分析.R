rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\富集分析") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DOSE")
#library(DOSE)

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("ggnewscale")
library("DOSE")
library(stringr)

pvalueFilter=0.05        
qvalueFilter=1        
showNum=20

rt=read.table("2.DIFF_all.txt",sep="\t",check.names=F,header=T)  
rownames(rt)=rt[,1]
rt=rt[read.table("disease.txt",sep="\t",check.names=F,header=F)[,1] ,c(1,2)]

genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","logFC","entrezID") 
#查看转换失败的ID
rt[rt[,"entrezID"]=="NA",]
rt=rt[rt[,"entrezID"]!="NA",]                        
gene=rt$entrezID
gene=unique(gene)
aflogfc=rt$logFC
names(aflogfc)=rt$symbol
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

pdf(file="KEGG_barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel) +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_cnet.pdf",width = 10,height = 8)
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(af, showCategory = 5, categorySize="pvalue",circular = TRUE,colorEdge = TRUE,cex_label_category=0.65,cex_label_gene=0.6,foldChange = aflogfc)
dev.off()

pdf(file="KEGG_net.pdf",width = 9,height = 7)
x2 <- pairwise_termsim(kk)
emapplot(x2,showCategory = showNum,cex_label_category=0.65,color = "pvalue",layout ="nicely")
dev.off()  


keggId="hsa05171"
geneFC=rt$logFC
names(geneFC)=gene
pv.out=pathview(gene.data = geneFC, pathway.id = keggId, species = "hsa", out.suffix = "pathview")
p <- pathview(gene.data = geneFC, pathway.id = keggId, species = "hsa", kegg.native = F, sign.pos="bottomleft", same.layer = F)



