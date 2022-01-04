setwd("~/work_space/4.ProjectET/batch31_rna")
library(ballgown)
library(pheatmap)
library(ggplot2)
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

#########  Rarg FPKM quantification; pData; DEG

data_directory = "/home1/gyang/work_space/4.ProjectET/batch31_rna/5.stringtie"
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
fwrite(gene_expression[,c(1:7,10,11,13)], "Rarg_KDKO_8Cexpr.csv",quote = F)



########################################
#################### quality control
########################################
#########PCA
gene_expression <- data.frame(fread("Rarg_KDKO_8Cexpr.csv"),check.names = F)[,1:7]
batch20 <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/batch20_rna/Rarg_Tfap2c_KD_expr.csv"),check.names = F)[,4:10]
batch15 <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/batch15_rna/Klf5_KD_expr.csv"),check.names = F)[,4:7]
normat <- merge(gene_expression,batch20,by="name")
normat <- merge(normat,batch15,by="name")
#two_c_gexpr <- two_c_gexpr[which(two_c_gexpr$name %in% normat$name),]
rownames(normat) <- normat$name
normat <- normat[,-1]
pca <-  prcomp(t(normat),center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pca_result<-as.data.frame(pca$x)
pca_result$condition <- c(rep("Ctrl",3),rep("Rarg_KDKO",3),rep("Rarg_KD",3),rep("Tfap2c_KD",3),rep("Klf5_KD",3))
pca_result$condition <- factor(pca_result$condition,levels = c("Ctrl","Rarg_KDKO","Rarg_KD","Tfap2c_KD","Klf5_KD"))
p <- ggplot(pca_result) +
  geom_point(aes(x=pca_result[,1],y=pca_result[,2],color=pca_result$condition),size=5)
p<-p+gran_theme+theme(legend.title =element_blank())+labs(x="PC1: 32 % variance",y="PC2: 24 % variance")
ggsave("Batch15_20_31_8Cembryo_KDKO_PCA.pdf",dpi=300)


########################################
#################### Repeats
########################################
batch31_repeats <- data.frame(fread("./3.2squire_count/8Cell_RargKDKO_subFcounts.txt"),check.names = F)[,c("Subfamily:Family:Class","fpkm" )]
colnames(batch31_repeats)[2]="8Cell_RargKDKO"
batch15_repeats <- data.frame(fread("~/work_space/4.ProjectET/batch15_rna/3.2squire_count/8Cell_Ctrl_subFcounts.txt"),check.names = F)[,c("Subfamily:Family:Class","fpkm")]
colnames(batch15_repeats)[2]="8Cell_Ctrl"
all_repeats <- merge(batch31_repeats,batch15_repeats,by="Subfamily:Family:Class")
all_repeats[,2:3] <- log2(all_repeats[,2:3]+1)
all_repeats$log2fc <- all_repeats$`8Cell_RargKDKO` - all_repeats$`8Cell_Ctrl`
