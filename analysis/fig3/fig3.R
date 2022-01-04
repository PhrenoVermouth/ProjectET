library(data.table)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(stringr)
library(tidyr)
setwd("~/work_space/4.ProjectET/analysis/fig3")

############
##### 1. PE TF: bubble plot
############
setwd("~/work_space/4.ProjectET/analysis/fig3/3.motif_bubble")
motif_gexpr <- data.frame(fread("motif_gexpr.txt"),check.names = F)
motif_gexpr$x1 = "a"

enh <- motif_gexpr[1:11,1:5]
enh$Gene <- factor(enh$Gene,levels = rev(enh$Gene))

p <- ggplot(enh,aes(x=x1, y=Gene)) +geom_point(aes(size=gexpr ,colour=pvalue)) 
p +  scale_color_gradient(low = "DarkSlateBlue", high = "Yellow3") + theme_bw()
ggsave("../Enhancer_cofactor_bubble.pdf",width = 3.5,height = 6)

pro <- motif_gexpr[12:21,1:5]
pro$Gene <- factor(pro$Gene,levels = rev(pro$Gene))

p <- ggplot(pro,aes(x=x1, y=Gene)) +geom_point(aes(size=gexpr ,colour=pvalue)) 
p +  scale_color_gradient(low = "DarkSlateBlue", high = "Yellow3") + theme_bw()
ggsave("../Promoter_cofactor_bubble.pdf",width = 3.5,height = 6)


############ selected_TF_expr
selected_motif_gexpr <- data.frame(fread("selected_motif_prolonged_gexpr.txt"),check.names = F)
colnames(selected_motif_gexpr) <- c("gene","D10.oocyte" ,"D14.oocyte" ,  "8w.oocyte" ,  "MII.oocyte"  , "PN5.zygote" ,  "Early.2.cell", "Late.2.cell" , "4.cell", "8.cell")
rownames(selected_motif_gexpr) <- selected_motif_gexpr$gene
selected_motif_gexpr <- selected_motif_gexpr[,-1]
p <- pheatmap(selected_motif_gexpr,cluster_rows = T,cluster_cols = F,color = colorRampPalette(colors = c("NavyBlue","black","Yellow1"))(100))
ggsave("../selected_TF_expr_heatmap.pdf",p$gtable,dpi=300,width=5,height=4)
# selected_motif_gexpr <- gather(selected_motif_gexpr,"stage","fpkm",-"gene")
# 
# p <- ggline(selected_motif_gexpr, "stage", "fpkm",
#             color = "gene", palette =  "paired")
# p<- p+ scale_color_npg()
# ggpar(p,x.text.angle = 90)
# ggsave("genomecov_aar_arr.pdf",dpi=300,width=3,height=4)





############
##### 2. PE-TT occupancy
############
setwd("~/work_space/4.ProjectET/analysis/fig3/4.PEloop_TToccu")
TToccu <- data.frame(fread("stat_on_TToccu_PE.txt"),check.names = F)
TToccu <- gather(TToccu,"type","num",-"class")
TToccu$type <- factor(TToccu$type,levels = c( "co_occu","only_Occu_by_Tfap2c","only_Occu_by_Tead4","none"))
ggbarplot(TToccu, "class", "num",
          fill = "type", color = "black", palette = c("#B22222","#FA8072","#FFC0CB","#FFFAFA"),
          label = TRUE, lab.col = "black", lab.pos = "in",width=0.5)
ggsave("../PEloop_TToccu_barplot.pdf",width = 3,height = 5)

