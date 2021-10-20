setwd("~/work_space/4.ProjectET/batch20_rna")
library(ballgown)
library(pheatmap)
library(ggplot2)
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

#########  Rarg FPKM quantification; pData; DEG

data_directory = "/home1/gyang/work_space/4.ProjectET/batch20_rna/5.stringtie/Rarg_KD"
bg = ballgown(dataDir=data_directory, samplePattern='rep', meas='all')
structure(bg)$exon
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(0,1), each=3))
phenotype_table = pData(bg)
gene_expression = data.frame(log2(gexpr(bg)+1),check.names = F)
gene_expression$name <- rownames(gene_expression)
colnames(gene_expression) <- gsub("FPKM.","",colnames(gene_expression))

gene_expression$Ctrl <- rowMeans(gene_expression[,1:3])
gene_expression$KD <- rowMeans(gene_expression[,4:6])
gene_expression$diff <- gene_expression$KD - gene_expression$Ctrl


stat_results = stattest(bg, feature="gene", meas='FPKM', covariate='group',libadjust=F)
stat_results <- stat_results[,c("id","pval","qval")]
#head(stat_results)
gene_expression <- merge(gene_expression, stat_results, by.x = "name", by.y="id")
gene_expression$class <- ifelse(abs(gene_expression$diff) > log2(1.5) & gene_expression$pval < 0.05,
                                ifelse(gene_expression$diff > 0, "up","down"),"none")
fwrite(gene_expression[,c(1:7,10,11,13)], "Rarg_KD_8Cexpr.csv",quote = F)

######### Tfap2c FPKM quantification; pData; DEG

data_directory = "/home1/gyang/work_space/4.ProjectET/batch20_rna/5.stringtie/Tfap2c_KD"
bg = ballgown(dataDir=data_directory, samplePattern='rep', meas='all')
structure(bg)$exon
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(0,1), each=3))
phenotype_table = pData(bg)
gene_expression = data.frame(log2(gexpr(bg)+1),check.names = F)
gene_expression$name <- rownames(gene_expression)
colnames(gene_expression) <- gsub("FPKM.","",colnames(gene_expression))

gene_expression$Ctrl <- rowMeans(gene_expression[,1:3])
gene_expression$KD <- rowMeans(gene_expression[,4:6])
gene_expression$diff <- gene_expression$KD - gene_expression$Ctrl


stat_results = stattest(bg, feature="gene", meas='FPKM', covariate='group',libadjust=F)
stat_results <- stat_results[,c("id","pval","qval")]
#head(stat_results)
gene_expression <- merge(gene_expression, stat_results, by.x = "name", by.y="id")
gene_expression$class <- ifelse(abs(gene_expression$diff) > log2(1.5) & gene_expression$pval < 0.05,
                                ifelse(gene_expression$diff > 0, "up","down"),"none")
fwrite(gene_expression[,c(1:7,10,11,13)], "Tfap2c_KD_8Cexpr.csv",quote = F)


########################################
#################### quality control
########################################
#########PCA
gene_expression <- data.frame(fread("Rarg_Tfap2c_KD_expr.csv"),check.names = F)
batch15 <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/batch15_rna/Klf5_KD_expr.csv"),check.names = F)
normat <- merge(gene_expression,batch15,by="name")
rownames(normat) <- normat$name
normat <- normat[,-1]
normat <- normat[,1:15]
pca <-  prcomp(t(normat),center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pca_result<-as.data.frame(pca$x)
pca_result$condition <- c(rep("Ctrl_batch20",3),rep("Rarg_KD",3),rep("Tfap2c_KD",3),rep("Ctrl_batch15",3),rep("Klf5_KD",3))
pca_result$condition <- factor(pca_result$condition,levels = c("Ctrl_batch20","Ctrl_batch15","Rarg_KD","Tfap2c_KD","Klf5_KD"))
p <- ggplot(pca_result) +
  geom_point(aes(x=pca_result[,1],y=pca_result[,2],color=pca_result$condition),size=5)
p<-p+gran_theme+theme(legend.title =element_blank())+labs(x="PC1: 37 % variance",y="PC2: 14 % variance")
ggsave("Batch15_20_8Cembryo_KD_PCA.pdf",dpi=300)
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


tfa <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/batch20_rna/Tfap2c_KD_8Cexpr.csv"),check.names = F)
tfa_up <- tfa[which(tfa$class=="up"),"name"] #748
tfa_down <- tfa[which(tfa$class=="down"),"name"] #869

klf <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/batch15_rna/FPKM-Klf5KD-DEG-by1.5-0.05.csv"),check.names = F)
klf_up <- klf[which(klf$class=="up"),"name"] #247
klf_down <- klf[which(klf$class=="down"),"name"] #574


intersect(tfa_up,klf_up) #66
intersect(tfa_down,klf_down) #187

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

twosetfisher(tfa_down,klf_down,24225,"Tfap2c_Klf5_kd_down.pdf","Tfap2c_Klf5_kd_down (< 2.2e-16)")


##### Tfap2c KD VS ICM TE
geneexp <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Normal_embryo_noem2c_log2FPKM.csv"),check.names = F)
boxplot(geneexp[which(geneexp$gene %in% tfa_down),c("ICM","TE")])
