
rm(list=ls())
options(stringsAsFactors = F)

setwd("D:\\肝细胞癌\\12.hub基因TCGA生存分析\\3.survival") 

#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)

inputFile="hub_survival.txt"      #输入文件
col=c("red","blue")          #定义高低表达组的颜色

rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt$futime=rt$futime/365    

#对基因进行循环
outTab=data.frame()
for(gene in colnames(rt)[3:ncol(rt)]){
  group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
  diff=survdiff(Surv(futime, fustat) ~group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outVector=cbind(gene,pValue)
  outTab=rbind(outTab,outVector)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     pval=pValue,
                     pval.size=6,
                     legend.labs=c("high","low"),
                     legend.title=paste0(gene," levels"),
                     font.legend=12,
                     xlab="futime(months)",
                     palette=col,
                     break.futime.by = 50,
                     conf.int=T,
                     fontsize=4.5,
                     risk.table=TRUE,
                     ylab="Overall survival",
                     risk.table.title="",
                     risk.table.height=.25)
  pdf(file=paste0(gene,".TCGA.pdf"),onefile = FALSE,
      width = 6,         #图片的宽度
      height =5)         #图片的高度
  print(surPlot)
  dev.off()
}
#输出基因和p值表格文件
write.table(outTab,file="TCGA.result.xls",sep="\t",row.names=F,quote=F)


