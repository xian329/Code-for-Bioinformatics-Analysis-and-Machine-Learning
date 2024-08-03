
rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\15.nomo+ROC-GEO") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DOSE")
#library(DOSE)


#引用包
library(dplyr)
library(pROC)
library(ggplot2)
library(survival)
library(regplot)
library(rms)
library(ggsci)
library(survminer)
library(timeROC)
library(ggDCA)
library(limma)
library(rmda)

inputFile="GSE84402.txt"       #表达矩阵
hub="hub.txt"        #核心基因

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#如果数据没有取log2,会对数据自动取log2
qx=as.numeric(quantile(data, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  data[data<0]=0
  data=log2(data+1)}
data2=normalizeBetweenArrays(data)
data2=data2-apply(data2,1,min)
data2=t(data2)

sample=read.table("sample.txt",sep="\t",header=F,check.names=F)
colnames(sample)=c("ID","Type")
data2=data2[sample$ID,]
HCC1=data2[,read.table(hub, header=F, sep="\t", check.names=F)[,1]]
HCC=cbind(sample,HCC1)
#简单看一下ROC曲线AUC的情况
aflist=roc(Type~SMG9+DDX49+NR2C2AP+BSG+OIT3+TPM2, data = HCC)
g3 <- ggroc(aflist, size = 1.2,alpha=.6,)
g5=g3+ggsci::scale_color_lancet()
print(g5)

#诺曼图，高度比：8比10
#Logistic回归模型
dd <- datadist(HCC)
options(datadist="dd")

fit <- lrm(formula = Type ~ SMG9+DDX49+NR2C2AP+BSG+OIT3+TPM2, 
           maxit = 1000, data =HCC,x=T, y=T)
print(fit)

coef=as.data.frame(fit$coefficients)[-1,,drop=F]
coefout=cbind(ID=rownames(coef),coef)
write.table(coefout,file="coefficients.txt",sep="\t",quote=F,row.names = F)

#绘图
pdf(file="nomogram.pdf", width=12, height=7.5)
#plot(nomogram(fit,fun.at = seq(0.05,0.95,0.05)),funlabel = "nomogram model")
nomo=nomogram(fit, fun=plogis,
              fun.at=c(0.0001,0.1,0.9,0.99),
              lp=F, funlabel="Risk of Disease")
plot(nomo)
dev.off()

#绘制校准曲线
cali=calibrate(fit, method="boot", B=1000)
pdf("hub_Calibration.pdf", width=5.5, height=5.5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()


nomoscore=predict(fit, data=t(HCC))
HCC$nomoscore=nomoscore
write.table(HCC,file="nomoscore.txt",sep="\t",quote=F,row.names = F)



