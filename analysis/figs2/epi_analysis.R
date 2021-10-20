setwd("~/work_space/4.ProjectET/analysis/figs3")
library(data.table)
library(stringr)

  
allele_df <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/project_et_allele.counts"),check.names = F)
allele_df <- allele_df[which(allele_df$chr != "chrX" ),]
allele_df <- allele_df[which(allele_df$chr != "chrY" ),]
allele_df <- allele_df[which(allele_df$chr != "chrM" ),]

temp <- allele_df[which(allele_df$`8Cell_Tead4_chr1-19_c57` + allele_df$`8Cell_Tead4_chr1-19_pwk` > 100),]
temp_c57 <- temp[which(temp$`8Cell_Tead4_chr1-19_c57` > 3*temp$`8Cell_Tead4_chr1-19_pwk` ),]
temp_pwk <- temp[which(temp$`8Cell_Tead4_chr1-19_pwk` > 3*temp$`8Cell_Tead4_chr1-19_c57` ),]
fwrite(temp_c57[,c("chr","start","end")],file = "8Cell_Tead4_c57.bed_temp", quote = F,sep="\t",col.names = F)
fwrite(temp_pwk[,c("chr","start","end")],file = "8Cell_Tead4_pwk.bed_temp", quote = F,sep="\t",col.names = F)

temp <- allele_df[which(allele_df$`8Cell_Tfap2c_chr1-19_c57` + allele_df$`8Cell_Tfap2c_chr1-19_pwk` > 100),]
temp_c57 <- temp[which(temp$`8Cell_Tfap2c_chr1-19_c57` > 3*temp$`8Cell_Tfap2c_chr1-19_pwk` ),]
temp_pwk <- temp[which(temp$`8Cell_Tfap2c_chr1-19_pwk` > 3*temp$`8Cell_Tfap2c_chr1-19_c57` ),]
fwrite(temp_c57[,c("chr","start","end")],file = "8Cell_Tfap2c_c57.bed_temp", quote = F,sep="\t",col.names = F)
fwrite(temp_pwk[,c("chr","start","end")],file = "8Cell_Tfap2c_pwk.bed_temp", quote = F,sep="\t",col.names = F)

temp <- allele_df[which(allele_df$`Morula_Tead4_chr1-19_c57` + allele_df$`Morula_Tead4_chr1-19_pwk` > 100),]
temp_c57 <- temp[which(temp$`Morula_Tead4_chr1-19_c57` > 3*temp$`Morula_Tead4_chr1-19_pwk` ),]
temp_pwk <- temp[which(temp$`Morula_Tead4_chr1-19_pwk` > 3*temp$`Morula_Tead4_chr1-19_c57` ),]
fwrite(temp_c57[,c("chr","start","end")],file = "Morula_Tead4_c57.bed_temp", quote = F,sep="\t",col.names = F)
fwrite(temp_pwk[,c("chr","start","end")],file = "Morula_Tead4_pwk.bed_temp", quote = F,sep="\t",col.names = F)

temp <- allele_df[which(allele_df$`Morula_Tfap2c_chr1-19_c57` + allele_df$`Morula_Tfap2c_chr1-19_pwk` > 100),]
temp_c57 <- temp[which(temp$`Morula_Tfap2c_chr1-19_c57` > 3*temp$`Morula_Tfap2c_chr1-19_pwk` ),]
temp_pwk <- temp[which(temp$`Morula_Tfap2c_chr1-19_pwk` > 3*temp$`Morula_Tfap2c_chr1-19_c57` ),]
fwrite(temp_c57[,c("chr","start","end")],file = "Morula_Tfap2c_c57.bed_temp", quote = F,sep="\t",col.names = F)
fwrite(temp_pwk[,c("chr","start","end")],file = "Morula_Tfap2c_pwk.bed_temp", quote = F,sep="\t",col.names = F) 
