
rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\疾病特征基因\\1.SVM-REF特征选择") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("sigFeature")
#library(DOSE)

#引用包
library(tidyverse)
library(glmnet)
source('msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
library(limma)

#输入文件
inputFile="symbol.txt"    # 表达矩阵   

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
#控制组放置最前
#分组
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
afcon=as.matrix(table(sample[,1]))[1,1]
afcon=as.vector(afcon)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(as.numeric(group))
rownames(group)=rownames(data)
colnames(group)="Type"
input <- as.data.frame(cbind(group,data))
input$Type=as.factor(input$Type)
#采用十折交叉验证
svmRFE(input, k = 10, halve.above = 100) #分割数据，分配随机数
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) #特征选择
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)
#把SVM-REF找到的特征保存到文件，AvgRank为按 10 次折叠的平均排名排序
write.csv(top.features,"feature_svm.csv")

# 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
# 选前40个变量进行SVM模型构建，然后导入已经运行好的结果
featsweep = lapply(1:40, FeatSweep.wrap, results, input) 

# 画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
pdf("svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误率
dev.off()

pdf("svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()

# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors) 

