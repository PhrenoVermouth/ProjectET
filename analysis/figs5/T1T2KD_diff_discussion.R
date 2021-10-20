library(data.table)
library(pheatmap)
library(stringr)
library(VennDiagram)


twosetfisher <- function(set1,set2,total_num,venn_name,venn_main=""){
  
  temp_int <- intersect(set1,set2)
  temp_merge <- unique(sort(c(set1,set2)))
  temp_mat <-
    matrix(c(length(temp_int),  length(setdiff(set1,temp_int)),length(setdiff(set2,temp_int)), total_num-length(temp_merge)),
           nrow = 2,
           dimnames = list(Guess = c("in", "out"),
                           Truth = c("in", "out")))
  print(fisher.test(temp_mat))
  
  pdf(venn_name)
  
  D1<-venn.diagram(list(set1=set1,set2=set2),filename=NULL,lwd=1,lty=1,col=c('red','blue'),
                   fill=c('red','blue'),cat.col=c('red','blue'),rotation.degree=180,main = venn_main)
  grid.draw(D1)
  dev.off()
  
}

######## 1. Down-regulated genes of T1T2 KD --- Overlap

setwd("/home1/gyang/work_space/4.ProjectET/analysis/figs5")
TTkd <- data.frame(fread("build/T1T2_KD_normat.csv"))
T1_down <- data.frame(fread("build/T1_KD_8Cexpr.csv"))
T1_down <- T1_down[which(T1_down$class == "down"),][,"name"]
T2_down <- data.frame(fread("build/T2_KD_8Cexpr.csv"))
T2_down <- T2_down[which(T2_down$class == "down"),][,"name"]

twosetfisher(T1_down,T2_down,24225,"Tfap2c_Tead4_KD_down_Venn.pdf","Tfap2c_Tead4_KD_down(< 2.2e-16)")

######## 2. Down-regulated genes of T1T2 KD --- heatmap

T1_down_expr <- TTkd[which(TTkd$name %in% T1_down),]
T2_down_expr <- TTkd[which(TTkd$name %in% T2_down),]
#TT_down_expr <- rbind(T1_down_expr, T2_down_expr)
#TT_down_expr <- TT_down_expr[,str_detect(colnames(TT_down_expr),"_A") ]
#TT_down_expr[,2:13] <- scale(TT_down_expr[,2:13])
#TT_down_expr <- t(scale(t(TT_down_expr)))
boxplot(T1_down_expr[,2:13],las=2)

boxplot(T2_down_expr[,2:13],las=2)

pheatmap(TT_down_expr,cluster_rows = F,cluster_cols = F)
