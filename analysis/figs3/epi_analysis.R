setwd("~/work_space/4.ProjectET/analysis/figs3")
library(data.table)
library(stringr)

############
##### 1. Output allele-biased peaks
############
  
allele_df <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/project_et_t2_allele.counts"),check.names = F)
allele_df <- allele_df[which(allele_df$chr != "chrX" ),]
allele_df <- allele_df[which(allele_df$chr != "chrY" ),]
allele_df <- allele_df[which(allele_df$chr != "chrM" ),]
#2463178 Kb

for (i in 1:6) {
  name=gsub(".c57","",colnames(allele_df)[2*i+2])
  name=paste0(name,".ASraw.txt")
  df <- allele_df[,c(1,2,3,2*i+2,2*i+3)]
  df$total <- df[,4] + df[,5]
  df$lfc <- abs(log2(df[,4]+1) - log2(df[,5]+1))
  df <- df[which(df$total >= 100),]
  fwrite(df[,1:7],name,quote = F,sep = "\t",col.names = F)
}

library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ChIPseeker)
txdb<- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/work_space/4.ProjectET/analysis/figs3/AS_stat")
peak_list <- list.files(pattern="peak.bed")

out_supp <- data.frame()
for (i in peak_list){
  name = gsub(".peak.bed","",i)
peakAnno <- data.frame(annotatePeak(i, tssRegion=c(-2500, 500),
                         TxDb=txdb, annoDb="org.Mm.eg.db"))
peakAnno <- peakAnno[,c("seqnames","start","end","SYMBOL", "annotation")]
peakAnno$origin = name
out_supp <- rbind(out_supp,peakAnno)
}
fwrite(out_supp,"../Allele_biased_peaks.txt",quote = F,sep = "\t",col.names = F)


############
##### 2. ABP versus SNP？
############
setwd("~/work_space/4.ProjectET/analysis/figs3/tf_snp/1.ABP_loci/check_snp_pm")

motif_snp_total <- data.frame(fread("tf_motif_snp_stat_tib.txt"),check.names = F)
motif_snp_total$xaxis <- paste0(motif_snp_total$tf,"_",motif_snp_total$allelic_type)
motif_snp_total$xaxis <- factor(motif_snp_total$xaxis,levels = c("Tfap2c_maternal","Tfap2c_paternal","Tead4_maternal","Tead4_paternal"))

p <- ggbarplot(motif_snp_total, "xaxis", "num",
               fill = "snp_type", color = "snp_type", palette = c("IndianRed1","SteelBlue"),width = 0.5)
p <- ggpar(p,legend = "none",x.text.angle = 45)
ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_summarybar.pdf",width = 3,height = 5.5)

motif_snp_logo <- data.frame(fread("TF_SNP_position_stat.txt"),check.names = F)
temp <- motif_snp_logo[motif_snp_logo$tf == "Tfap2c",]
# p <- ggline(temp, "position", "snp_percent",
#             color = "allelic_type", palette =  c("IndianRed1","SteelBlue"))
# ggpar(p,width=6,height=4)
# ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_line_tfap2c.pdf",width = 5,height = 3)

temp$allelic_type <- factor(temp$allelic_type,levels = c("maternal","paternal"))
p <- ggbarplot(temp, "position", "snp_percent",
               fill = "allelic_type", color = "white", palette = c("#282A73","#96CB6B"),alpha=0.8,width = 0.5,position = position_dodge(0.5))
p <- ggpar(p,legend = "none")
ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_bar_tfap2c.pdf",width = 5,height = 3)

motif_snp_logo <- data.frame(fread("TF_SNP_position_stat.txt"),check.names = F)
temp <- motif_snp_logo[motif_snp_logo$tf == "Tead4",]
# p <- ggline(temp, "position", "snp_percent",
#             color = "allelic_type", palette =  c("IndianRed1","SteelBlue"))
# ggpar(p,width=6,height=4)
# ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_line_tead4.pdf",width = 5,height = 3)
temp$allelic_type <- factor(temp$allelic_type,levels = c("maternal","paternal"))
p <- ggbarplot(temp, "position", "snp_percent",
               fill = "allelic_type", color = "white", palette = c("#282A73","#96CB6B"),alpha=0.8,width = 0.5,position = position_dodge(0.5))
p <- ggpar(p,legend = "none")
ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_bar_tead4.pdf",width = 5,height = 3)