
rm(list=ls())#清空环境以排除干扰
options(stringsAsFactors = F)#防止文件乱码

library(dplyr)
library(tidyverse)

setwd("D:\\肝细胞癌\\11.TCGA临床信息提取")

#从clinical中筛选临床特征
clinical <- read_tsv("clinical.tsv")
colnames(clinical)
clin_data <- clinical %>% 
  dplyr::select(patient_ID = case_submitter_id,
                vital_status,
                age = age_at_diagnosis,
                gender = gender,
                Stage = ajcc_pathologic_stage,
                days_to_death,
                days_to_last_follow_up) 

clin_data$age <- (as.numeric(clin_data$age))/365
clin_data$days_to_death <- as.numeric(clin_data$days_to_death)
clin_data$days_to_last_follow_up <- as.numeric(clin_data$days_to_last_follow_up)

clin_data1 <- clin_data %>%
  mutate(futime = coalesce(days_to_death, days_to_last_follow_up),
         fustat = ifelse(vital_status == 'Alive', 0, 1)) %>%
  dplyr::filter(futime > 0.08) %>%
  distinct() %>%
  dplyr::select(- c("days_to_death", "days_to_last_follow_up")) %>% 
  na.omit()

write.table(clin_data1,"preclinical.txt",sep = "\t",row.names=FALSE,col.names=TRUE)


library(dplyr)
library(tidyverse)

rt1=read.table('preclinical.txt',sep="\t",header=T,check.names=F)

#Stage
rt1$Stage<-str_replace_all(rt1$Stage, c("Stage IIIA" = "Stage III"))
rt1$Stage<-str_replace_all(rt1$Stage, c("Stage IIIB" = "Stage III"))
rt1$Stage<-str_replace_all(rt1$Stage, c("Stage IIIC" = "Stage III"))
rt1$Stage<-str_replace_all(rt1$Stage, c("Stage IVB" = "Stage IV"))
rt1$Stage<-str_replace_all(rt1$Stage, c("Stage IVA" = "Stage IV"))
rt1$Stage<-str_replace_all(rt1$Stage, c("'--" = "unknow"))


#整理vital_status状况
table(rt1$vital_status)
rt1$vital_status <- str_replace_all(rt1$vital_status, c("Not Reported" = "unknow"))


#保存
write.table(rt1,"clinical_data.txt",sep = "\t",row.names = F,quote = F)




