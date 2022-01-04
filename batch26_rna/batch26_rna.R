setwd("~/work_space/4.ProjectET/batch26_rna/5.stringtie")
library(ballgown)
library(pheatmap)
library(ggplot2)
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

#########  Rarg FPKM quantification; pData; DEG

data_directory = "~/work_space/4.ProjectET/batch26_rna/5.stringtie"
bg = ballgown(dataDir=data_directory, samplePattern='rep', meas='all')
structure(bg)$exon
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(0,1), each=3))
phenotype_table = pData(bg)
gene_expression = data.frame(log2(gexpr(bg)+1),check.names = F)
gene_expression$name <- rownames(gene_expression)
colnames(gene_expression) <- gsub("FPKM.","",colnames(gene_expression))

gene_expression$epi <- rowMeans(gene_expression[,1:3])
gene_expression$exe <- rowMeans(gene_expression[,4:5])
gene_expression$diff <- gene_expression$epi - gene_expression$exe


# fwrite(gene_expression[,c(1:6)], "Embryo_E65_RNA_yizhang_xw.csv",quote = F)
# fwrite(gene_expression[,c(6:8)], "Embryo_E65_RNA_yizhang_xw_collapse.csv",quote = F)
