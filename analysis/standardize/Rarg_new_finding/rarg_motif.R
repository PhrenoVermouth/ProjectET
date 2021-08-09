library(clusterProfiler)
library(ChIPseeker)
library(data.table)
library(ggpubr)
library(org.Mm.eg.db)
library(tidyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
setwd("~/work_space/4.ProjectET/analysis/standardize/Rarg_new_finding")

######## 5mC dynamics during ZGA

l1 <- data.frame(fread("build/l1"))
l1$zygote <- (l1$GSM1386019_oocyte + l1$GSM1386020_sperm)/2 
l1$tag <- "l1"

l2 <- data.frame(fread("build/l2"))
l2$zygote <- (l2$GSM1386019_oocyte + l2$GSM1386020_sperm)/2 
l2$tag <- "l2"

l3 <- data.frame(fread("build/l3"))
l3$zygote <- (l3$GSM1386019_oocyte + l3$GSM1386020_sperm)/2
l3$tag <- "l3"

l123 <- rbind(l1,l2)
l123 <- rbind(l123,l3)
colnames(l123) <- c("chr" , "start" ,"end" ,"oocyte","sperm","2cell", "4cell","zygote", "tag") 
#l123 <- l123[,c("chr" , "start" ,"end" ,"2cell", "4cell","zygote", "tag")]
l123 <- gather(l123,"time","ratio",-"chr",-"start",-"end",-"tag")
l123$time <- factor(l123$time,levels = c("oocyte","sperm","zygote","2cell", "4cell"))

p <- ggboxplot(l123,"time","ratio",facet.by = "tag")
ggpar(p,x.text.angle = 45)
ggsave("SINE_Rarg_3levels.png",dpi = 300,width = 7, height = 3)

######## l1 promoter?
#peak <- readPeakFile()
peakAnno <- annotatePeak("~/work_space/4.ProjectET/analysis/standardize/Rarg_new_finding/rarg_l1.bed", tssRegion=c(-2500, 500),
             TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_df <- data.frame(peakAnno)
unique(peakAnno_df$SYMBOL)
peakAnno_df <- peakAnno_df[which(peakAnno_df$annotation == "Promoter (<=1kb)"),]
sort(unique(peakAnno_df$SYMBOL))
fwrite(list(sort(unique(peakAnno_df$SYMBOL))),file = "L1_RArg_promoter1kb", quote = F,sep=",",col.names = T)
geneexp <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Xiewei_2016_Nature_log2FPKM_RNA.csv"))
geneexp[,2:12] <- scale(geneexp[,2:12]) 
boxplot(geneexp[which(geneexp$Gene %in% sort(unique(peakAnno_df$SYMBOL))),2:12],las=2)
