setwd("~/work_space/4.ProjectET/batch15_rna")
library(ballgown)
library(pheatmap)
library(ggplot2)
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ChIPseeker)
library(ggplot2)
library(org.Mm.eg.db)
txdb<- TxDb.Mmusculus.UCSC.mm10.knownGene

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


stat_results = stattest(bg, feature='gene', meas='FPKM', covariate='group',libadjust=F)
stat_results <- stat_results[,c("id","pval","qval")]
#head(stat_results)
gene_expression <- merge(gene_expression, stat_results, by.x = "name", by.y="id")
gene_expression$class <- ifelse(abs(gene_expression$diff) > log2(1.5) & gene_expression$pval < 0.05,
                                ifelse(gene_expression$diff > 0, "up","down"),"none")
fwrite(gene_expression[,c(1:7,10,13)], "Klf5_KD_220103altered_expr.csv",quote = F)



key_influ <- c("Arpc1b","Zscan4c","Zscan4d","Zscan4f","Zfp532","Zfp560","Dppa3")
all_fpkm2$highlight <- ifelse(all_fpkm2$gene %in% key_influ, "yes", "no")
all_fpkm2 <- spread(all_fpkm2,"class","fpkm")
p <- ggscatter(all_fpkm2, x = "IVF2C", y = "SCNT2C",alpha = 1,color = "highlight", palette = c("grey","Firebrick4"))+
  geom_abline(linetype = "dashed",intercept = 1,alpha=0.8,color = "black") +
  geom_abline(linetype = "dashed",intercept = -1,alpha=0.8,color = "black") +
  xlim(-1,18) + ylim(-1,18)+coord_cartesian(expand = F)




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

#########  Scatter plot

key_influ <- c("Actn1","Actn4","Pard6b","Plcd1","Plcg1")
all_fpkm2 <- gene_expression[,c("name","Ctrl","KD","diff","pval", "class")]
all_fpkm2$highlight <- ifelse(all_fpkm2$name %in% key_influ, "yes", "no")
#all_fpkm2 <- spread(all_fpkm2,"class","fpkm")
p <- ggscatter(all_fpkm2, x = "Ctrl", y = "KD",alpha = 1,color = "class", palette = c("#4682B4","grey","#FF6A6A"))+
  geom_abline(linetype = "dashed",intercept = 0.585,alpha=0.8,color = "black") +
  geom_abline(linetype = "dashed",intercept = -0.585,alpha=0.8,color = "black") +
  xlim(0,12) + ylim(0,12)+coord_cartesian(expand = F)

ggsave("~/work_space/4.ProjectET/analysis/fig3/5.Klf_enhancer_binding/Klf5KD_scatter_DEG.pdf",width = 4,height = 4.5,useDingbats=F) 

gene_expression[gene_expression$class == "down","name"]

ego2 <- enrichGO(gene         = gene_expression[gene_expression$class == "down","name"],
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dotplot(ego2, showCategory=10) + ggtitle("GO for Klf5KD Down-regulated genes") 

ggsave("~/work_space/4.ProjectET/analysis/fig3/Klf5KD_DEG_downGO.pdf",width=7 ,height = 4)

