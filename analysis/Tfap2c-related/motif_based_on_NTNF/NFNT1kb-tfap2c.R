raw_density <- data.frame(fread("~/work_space/1.Mouse_Acetylation/buffet/k9_analysis/2.merged_bw/unjusted_bw/morula_nf_nt_1kb.tab"),check.names = F)
raw_density[,c(4:ncol(raw_density))] <- log2(raw_density[,c(4:ncol(raw_density))] +1) 

morula_density <- raw_density[,c("chr","start","end","K9-morula","K9-NT-morula")]
morula_density$fc <- morula_density$`K9-NT-morula` - morula_density$`K9-morula`

quantile(morula_density$fc,probs = seq(0, 1, 0.02)) #-0.97 2%
morula_density <- morula_density[which(morula_density$fc < -0.97),]
morula_density <- morula_density[which(morula_density$chr != "chrX"),]
morula_density <- morula_density[which(morula_density$chr != "chrY"),]
fwrite(morula_density[,c(1,2,3)],file = "~/work_space/1.Mouse_Acetylation/buffet/k9_analysis/2.merged_bw/unjusted_bw/test.bed", quote = F,sep="\t",col.names = F)


morula_density$color <- ifelse(morula_density$fc >= 1.5,"NT_specific",ifelse(morula_density$fc <= -1.5,"WT_specific","plain"))
morula_density$color <- factor(morula_density$color ,levels=c("NT_specific","plain","WT_specific"))


fwrite(wt_nt_expr[which(wt_nt_expr$gene %in% c("Tfap2c","Tfap2a","Gsc","Tead3","Elf5","Tead4","Tead1","Crx","Gata3","Elf4","Pitx1","Klf5")),],file = "~/work_space/1.Mouse_Acetylation/buffet/k9_analysis/2.merged_bw/unjusted_bw/test.csv", quote = F,sep=",",col.names = T)

temp <- wt_nt_expr[which(wt_nt_expr$gene %in% c("Vangl1", "Vav1", "Trip10", "Tpm4", "Tmsb10", "Tbcd", "Tagln2", "Plekhg2", "Nck2", "Mylpf", "Marcksl1", "Marcks", "Lcp1", "Frmpd1", "Frmd4b", "Fgfr2", "Epb4.1l5", "Cdc42ep4", "Cdc42ep3", "Cdc42ep1", "Cdc42ep1", "Arpc1b", "Arhgef19", "Arhgef16", "Arhgdib", "Arhgap9", "Arhgap18", "Anln", "Amn")),c("gene","Morula","NT-morula")]
rownames(temp) <- temp$gene
temp <- temp[c("Vangl1", "Vav1", "Trip10", "Tpm4", "Tmsb10", "Tbcd", "Tagln2", "Plekhg2", "Nck2", "Mylpf", "Marcksl1", "Marcks", "Lcp1", "Frmpd1", "Frmd4b", "Fgfr2", "Epb4.1l5", "Cdc42ep4", "Cdc42ep3", "Cdc42ep1", "Cdc42ep1", "Arpc1b", "Arhgef19", "Arhgef16", "Arhgdib", "Arhgap9", "Arhgap18", "Anln", "Amn"),]
temp <- temp[,c(2,3)]
pheatmap::pheatmap(temp,cluster_rows = F)


#################### Tfap2c morula
### Homer motif locus calc

setwd("~/work_space/4.ProjectET/analysis/Tfap2c-related/tf_ana")
raw_bed <- fread("~/work_space/4.ProjectET/analysis/Tfap2c-related/test.bed")
raw_bed$pos <- as.numeric(rownames(raw_bed))
rela_locus <- fread("Tfap2c_1k.peak")
bed <- merge(raw_bed,rela_locus, by.x = c("pos"), by.y = c("PositionID"))
bed$motif_fake_start <- bed$V2 + bed$Offset + 1 #fake means suitable for + strand
bed$motif_fake_end <- bed$motif_fake_start + nchar(bed$Sequence) -1
bed$start <- ifelse(bed$Strand == "+", bed$motif_fake_start, 2*bed$motif_fake_start - bed$motif_fake_end)
bed$end <- ifelse(bed$Strand == "+", bed$motif_fake_end, bed$motif_fake_start)
# bed <- merge(bed, promoter, by.x = c("V1","V2"), by.y = c("V1","V2"))[,c("V1","motif_start","motif_end","Sequence","V4")]
# colnames(bed) = c("chr","start","end","seq","gene")
write.table(bed[,c("V1","start","end")],"Tfap2c_NFNT1kb.bed",sep="\t",quote = F,row.names = F,col.names = F)
