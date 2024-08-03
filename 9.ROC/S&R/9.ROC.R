
rm(list=ls())#清空环境以排除干扰
options(stringsAsFactors = F)#防止文件乱码

##ROC分析目标基因在TCGA_LIHC中的预测价值

#BiocManager::install("pROC")

library(pROC)

setwd("D:\\肝细胞癌\\9.ROC")  

rt=read.table("ROC_input.txt",header=T,sep="\t",check.names=F)

p1=roc(Type ~ TPM2, rt)
p2=roc(Type ~ DMGDH, rt)
p3=roc(Type ~ ACADSB, rt)
p4=roc(Type ~ OIT3, rt)
p5=roc(Type ~ GBA3, rt)
p6=roc(Type ~ nomoscore, rt)
#p7=roc(Type ~ AAAA, rt)

#绘制ROC曲线，col控制roc曲线的颜色
p1$auc
p2$auc
p3$auc
p4$auc
p5$auc
p6$auc
#p7$auc

pdf(file="TPM2.pdf",width=5,height=5)
plot(p1,
     print.auc=TRUE,
     auc.polygon=TRUE,
     grid=c(0.1,0.2),
     grid.col=c("green", "red"),
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue",
     print.thres=TRUE)
dev.off()

pdf(file="DMGDH.pdf",width=5,height=5)
plot(p2,
     print.auc=TRUE,
     auc.polygon=TRUE,
     grid=c(0.1,0.2),
     grid.col=c("green", "red"),
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue",
     print.thres=TRUE)
dev.off()

pdf(file="ACADSB.pdf",width=5,height=5)
plot(p3,
     print.auc=TRUE,
     auc.polygon=TRUE,
     grid=c(0.1,0.2),
     grid.col=c("green", "red"),
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue",
     print.thres=TRUE)
dev.off()

pdf(file="OIT3.pdf",width=5,height=5)
plot(p4,
     print.auc=TRUE,
     auc.polygon=TRUE,
     grid=c(0.1,0.2),
     grid.col=c("green", "red"),
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue",
     print.thres=TRUE)
dev.off()

pdf(file="GBA3.pdf",width=5,height=5)
plot(p5,
     print.auc=TRUE,
     auc.polygon=TRUE,
     grid=c(0.1,0.2),
     grid.col=c("green", "red"),
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue",
     print.thres=TRUE)
dev.off()


pdf(file="nomoscore.pdf",width=5,height=5)
plot(p6,
     print.auc=TRUE,
     auc.polygon=TRUE,
     grid=c(0.1,0.2),
     grid.col=c("green", "red"),
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue",
     print.thres=TRUE)
dev.off()

#pdf(file="AAAA.pdf",width=5,height=5)
#plot(p7,
#     print.auc=TRUE,
#     auc.polygon=TRUE,
#     grid=c(0.1,0.2),
#     grid.col=c("green", "red"),
#     max.auc.polygon=TRUE, 
#     auc.polygon.col="skyblue",
#     print.thres=TRUE)
#dev.off()
