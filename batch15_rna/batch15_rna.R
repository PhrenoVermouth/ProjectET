setwd("~/work_space/4.ProjectET/batch15_rna")
library(ballgown)
library(pheatmap)
library(ggplot2)
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

#########  FPKM quantification

# data_directory = "/home1/gyang/work_space/4.ProjectET/batch15_rna/5.stringtie"
# bg = ballgown(dataDir=data_directory, samplePattern='rep', meas='all')
# structure(bg)$exon
# gene_expression = data.frame(log2(gexpr(bg)+1),check.names = F)
# gene_expression$name <- rownames(gene_expression)
# colnames(gene_expression) <- gsub("FPKM.","",colnames(gene_expression))
# 
# gene_expression$mean1 <- rowMeans(gene_expression[,1:3])
# gene_expression$mean2 <- rowMeans(gene_expression[,4:6])
# gene_expression$diff <- gene_expression$mean1 - gene_expression$mean2
# 
# fwrite(gene_expression[,c(1:7,10)], "Klf5_KD_expr.csv",quote = F)

data_directory = "/home1/gyang/work_space/4.ProjectET/batch15_rna/5.stringtie"
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


stat_results = stattest(bg, feature='gene', meas='FPKM', covariate='group')
stat_results <- stat_results[,c("id","pval","qval")]
#head(stat_results)
gene_expression <- merge(gene_expression, stat_results, by.x = "name", by.y="id")
gene_expression$class <- ifelse(abs(gene_expression$diff) > log2(1.5) & gene_expression$pval < 0.05,
                                ifelse(gene_expression$diff > 0, "up","down"),"none")
fwrite(gene_expression[,c(1:7,10,13)], "Klf5_KD_0914altered_expr.csv",quote = F)

#########  DEG
# 
# setwd("~/work_space/4.ProjectET/batch15_rna")
# library(DESeq2)
# rawmat <- read.csv("5.featurecounts/merge_geneCounts_featureCounts_hisat2.tab",sep="\t",header=T,check.names = F)
# rawmat <- rawmat[!duplicated(rawmat[,1]),]
# rownames(rawmat) <- rawmat$`#id`
# rawmat <- rawmat[,-1]
# rawmat <- rawmat[,-1]
# 
# coldata <- read.csv("coldata.csv",row.names = 1,header = T,strip.white = T)
# all(rownames(coldata) %in% colnames(rawmat))
# all(rownames(coldata) == colnames(rawmat))
# rawmat <- rawmat[,rownames(coldata)]
# all(rownames(coldata) == colnames(rawmat))
# 
# ####################Deseq
# dds <- DESeqDataSetFromMatrix(countData = rawmat,colData = coldata,design = ~ condition)
# keep <- rowSums(counts(dds)) >= 5
# dds <- dds[keep,]
# #dds$condition <- relevel(dds$condition, ref = "untreated")
# dds <- DESeq(dds)
# res <- results(dds) #FDR=0.1 by default
# resOrdered <- res[order(res$padj),]
# sum(abs(res$log2FoldChange) > 1, na.rm = TRUE)
# rlog <- rlog(dds,blind=F)
# #assay(rlog["2310015A10Rik"])
# ##################Normalized matrix
# normat <- counts(dds,normalized=T)
# normat <- log2(counts(dds,normalized=T)+1)
# normat <- data.frame(normat,check.names = F)
# #normat <- scale(normat)
# #write.csv(normat,"all_normat.csv",quote=F)
# 


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
p<-p+gran_theme+theme(legend.title =element_blank())+labs(x="PC1: 43 % variance",y="PC2: 16 % variance")
ggsave("Klf5_KD_0914altered_PCA.pdf",dpi=300)

### TT target gene on Science Fig5G
diff_fpkm_DEG[which(diff_fpkm_DEG$name %in% c("Vangl1", "Vav1", "Trip10", "Tpm4", "Tmsb10", "Tbcd", "Tagln2", "Plekhg2", "Nck2", "Mylpf", "Marcksl1", "Marcks", "Lcp1", "Frmpd1", "Frmd4b", "Fgfr2", "Epb4.115", "Cdc42ep4", "Cdc42ep3", "Cdc42ep1", "Cdc42ep1", "Arpc1b", "Arhgef19", "Arhgef16", "Arhgdib", "Arhgap9", "Arhgap18", "Anln", "Amn")),]

fwrite(diff_fpkm_DEG,"FPKM-Klf5KD-DEG-by1.5-0.05.csv")
