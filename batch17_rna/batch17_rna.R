setwd("~/work_space/4.ProjectET/batch17_rna")
library(ballgown)
library(pheatmap)
library(ggplot2)
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

#########  FPKM quantification; pData; DEG

data_directory = "/home1/gyang/work_space/4.ProjectET/batch17_rna/5.stringtie"
bg = ballgown(dataDir=data_directory, samplePattern='rep', meas='all')
structure(bg)$exon
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(0,1), each=2))
phenotype_table = pData(bg)
gene_expression = data.frame(log2(gexpr(bg)+1),check.names = F)
gene_expression$name <- rownames(gene_expression)
colnames(gene_expression) <- gsub("FPKM.","",colnames(gene_expression))

gene_expression$Ctrl <- rowMeans(gene_expression[,1:2])
gene_expression$KD <- rowMeans(gene_expression[,3:4])
gene_expression$diff <- gene_expression$KD - gene_expression$Ctrl


stat_results = stattest(bg, feature="gene", meas='FPKM', covariate='group',libadjust=F)
stat_results <- stat_results[,c("id","pval","qval")]
#head(stat_results)
gene_expression <- merge(gene_expression, stat_results, by.x = "name", by.y="id")
gene_expression$class <- ifelse(abs(gene_expression$diff) > log2(1.5) & gene_expression$pval < 0.05,
                                ifelse(gene_expression$diff > 0, "up","down"),"none")
fwrite(gene_expression[,c(1:5,8,9,11)], "Tfap2c_KD_expr_trimmed.csv",quote = F)


########################################
#################### quality control
########################################
#########PCA
normat <- gene_expression[,c(1:7)]
rownames(normat) <- normat$name
normat <- normat[,-1]
pca <-  prcomp(t(normat),center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pca_result<-as.data.frame(pca$x)
pca_result$condition <- c(rep("Ctrl",3),rep("KD",3))
pca_result$condition <- factor(pca_result$condition,levels = c("Ctrl","KD"))
p <- ggplot(pca_result) +
  geom_point(aes(x=pca_result[,1],y=pca_result[,2],color=pca_result$condition),size=5)
p<-p+gran_theme+theme(legend.title =element_blank())+labs(x="PC1: 72 % variance",y="PC2: 13 % variance")
ggsave("Tfap2c_KD_PCA.pdf",dpi=300)
# p + scale_color_manual(values=c("LightGray","Cyan","Blue","magenta","red"))
# ggsave("newdata_all_PCA_newcol.pdf",dpi=300)

# temp = data.frame(nt=apply(normat[,1:5],1,mean),tsa=apply(normat[,6:10],1,mean),row.names = rownames(normat))
# temp$fc = temp$tsa - temp$nt
# temp$gene = rownames(temp)
#
# temp2=data.frame(res)

########################################
#################### DEG expr in pre embryos
########################################

##### Down ---- expr in ICM and TE lineage

gene_expression[which(gene_expression$class == "down"),c("name")]
geneexp <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Normal_embryo_noem2c_log2FPKM.csv"),check.names = F)
#geneexp <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/XieWei_all_embryos_log2.csv"),check.names = F)
geneexp <- geneexp[,c(1:9)]
#geneexp <- geneexp[,c("Gene","MII.oocyte","PN5.zygote","Early.2.cell" ,"Late.2.cell" ,"4.cell","8.cell" ,"E35ICM" , "E35TE" )]
geneexp_deg <- geneexp[which(geneexp$gene %in% gene_expression[which(gene_expression$class == "down"),c("name")]),]
#geneexp_deg <- geneexp[which(geneexp$Gene %in% gene_expression[which(gene_expression$class == "down"),c("name")]),]
rownames(geneexp_deg) <- geneexp_deg$gene
#rownames(geneexp_deg) <- geneexp_deg$Gene
geneexp_deg <- geneexp_deg[,-1]
#geneexp_deg <- geneexp_deg[,c("8cell","Morula","ICM","TE")]
geneexp_deg$fc <-  geneexp_deg$TE - geneexp_deg$ICM
#geneexp_deg$fc <-  geneexp_deg$E35TE - geneexp_deg$E35ICM
geneexp_deg <- data.frame(geneexp_deg[order(geneexp_deg$fc,decreasing = T),],check.names = F)
fwrite(geneexp_deg, "Tfap2c_KD_downregulated_trimmed.csv",quote = F,row.names = T)
#fwrite(geneexp_deg, "Tfap2c_KD_downregulated_Xw.csv",quote = F,row.names = T)
geneexp_deg <- geneexp_deg[,c(1:8)]
geneexp_deg <- na.omit(t(scale(t(geneexp_deg))))

p <- pheatmap(geneexp_deg,cluster_cols = F,cluster_rows = F, show_rownames = F)
ggsave("siTfap2c_downregulate_trimmed.pdf",p$gtable,dpi=300)

##### Up ---- 2C



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

Tfap2c_kd_up <- gene_expression[which(gene_expression$class == "up"),c("name")]
two_c <- data.frame(fread("~/work_space/1.Mouse_Acetylation/buffet/resource_set/2cell_3stages_specific_genes.txt",check.names = F,header = F))[,1]

twosetfisher(Tfap2c_kd_up,two_c,24225,"Tfap2c_kd_up_2C.pdf","Tfap2c_kd_up_2C (< 2.2e-16)")

