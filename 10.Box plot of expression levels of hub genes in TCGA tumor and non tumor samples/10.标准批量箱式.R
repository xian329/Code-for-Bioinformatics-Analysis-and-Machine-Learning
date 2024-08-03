rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\10.hub基因在TCGA中表达差异箱式图") 

#引用包
library(ggpubr)          
library(limma)
library(ggsci)

#输入文件  
expfile="symbol.txt"          #表达矩阵
hub="hub.txt"                     #核心基因
sample="sample.txt"               #样本属性
Ccol="blue"                       #正常控制组颜色
Pcol="red"                        #疾病实验组颜色
afmethod="wilcox.test"        #t检验："t.test"（小样本量）,秩和检验:"wilcox.test"（大样本量）

rt=read.table(expfile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)

rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]

dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

#如果数据没有取log2,会对数据自动取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)
data=data-apply(data,1,min)

data=t(data)
data2=read.table(sample,sep="\t",header=F,check.names=F,row.names = 1)
ssample=as.character(intersect(rownames(data2),rownames(data)))
data=data[ssample,]
data2=data2[ssample,,drop=F]
data1=cbind(data2,data)
colnames(data1)[1]="Type"
data1$ID=rownames(data1)
genes=read.table(hub,sep="\t",header=F,check.names=F)[,1]
data1=data1[,c("ID","Type",genes)]
af=colnames(data1)
df=c(3:ncol(data1))

for (i in df) {

outFile=paste0(af[i],"_","boxplotdiff.pdf")             
data=as.data.frame(data1[,c(1,2,i)])
x=colnames(data)[2]
y=colnames(data)[3]
colnames(data)=c("id","Type","Expression")


group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=group)
data$Expression=as.numeric(data$Expression)

comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


boxplot=ggboxplot(data, x="Type", y="Expression", fill="Type",
		          xlab=x,     
		          ylab=y,
		          legend.title=x,
		          palette = c(Ccol,Pcol),   
		          #add = "jitter"
		          )+ 
	    stat_compare_means(comparisons = my_comparisons,method=afmethod)


pdf(file=outFile,width=4,height=4.5)
print(boxplot)
dev.off()

}


