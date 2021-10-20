library(data.table)
library(stringr)

####################################
###### STAT on TFs occupancy on genes
####################################

setwd("/home1/gyang/work_space/4.ProjectET/analysis/fig1/1.dynamic_change/1.peak_on_gene")
all_peak_gene <- list.files(pattern = "?gene_loc.txt$")

# stage_tf_class <- list(Eight_cell = list(Tfap2c=list(promoter=c(),distal=c()), Tead4=list(promoter=c(),distal=c())),
#                        Morula = list(Tfap2c=list(promoter=c(),distal=c()), Tead4=list(promoter=c(),distal=c())))
stage_tf_class <- list()
for (i in all_peak_gene){
  #i=all_peak_gene[1]
  i_split <- str_split(i,"_",5)[[1]]
  i_split=gsub("8Cell","Eight_cell",i_split)
  stage = i_split[1]
  tf = i_split[2]
  print(stage)
  print(tf)
  df_tmp <- data.frame(fread(i))
  promoter <- unique(df_tmp[which(df_tmp$V8 == "promoter"),"V4"])
  genebody <-  unique(df_tmp[which(df_tmp$V8 == "genebody"),"V4"])
  stage_tf_class[[stage]][[tf]][["promoter"]] <- promoter
  stage_tf_class[[stage]][[tf]][["genebody"]] <- genebody
}

length(intersect(stage_tf_class[["Eight_cell"]][["Tead4"]][["promoter"]],stage_tf_class[["Eight_cell"]][["Tead4"]][["genebody"]]))
length(intersect(stage_tf_class[["Eight_cell"]][["Tfap2c"]][["promoter"]],stage_tf_class[["Eight_cell"]][["Tfap2c"]][["genebody"]]))

length(intersect(stage_tf_class[["Eight_cell"]][["Tead4"]][["promoter"]],stage_tf_class[["Eight_cell"]][["Tfap2c"]][["promoter"]]))
length(setdiff(stage_tf_class[["Eight_cell"]][["Tead4"]][["promoter"]],stage_tf_class[["Eight_cell"]][["Tfap2c"]][["promoter"]]))
length(setdiff(stage_tf_class[["Eight_cell"]][["Tfap2c"]][["promoter"]],stage_tf_class[["Eight_cell"]][["Tead4"]][["promoter"]]))

############ Output for epi-profiling
######8c
tf_promoter <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/standardize/3.tf_promoter/tf_promoter.tab"),check.names = F)
tf_promoter[,4:ncol(tf_promoter)] <- log2(tf_promoter[,4:ncol(tf_promoter)] + 1)
promoter_anno <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/1.dynamic_change/mm10_2k_0.5k_promoter.bed"))
tf_promoter <- merge(tf_promoter,promoter_anno, by.x=c("chr","start","end"), by.y=c("V1","V2","V3"))
t1_dom <- tf_promoter[which(tf_promoter$`8Cell_Tfap2c` - tf_promoter$`8Cell_Tead4` > 1),"V4"]
t2_dom <- tf_promoter[which(tf_promoter$`8Cell_Tead4` - tf_promoter$`8Cell_Tfap2c` > 1),"V4"]

promoter_anno <-  data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/fig1/1.dynamic_change/mm10_2k_0.5k_promoter.bed"))

t1 <- stage_tf_class[["Eight_cell"]][["Tfap2c"]][["promoter"]]
t2 <- stage_tf_class[["Eight_cell"]][["Tead4"]][["promoter"]]
t1_and_t2 <- intersect(t1,t2)
t1_only <- intersect(setdiff(t1,t2),t1_dom)
t2_only <- intersect(setdiff(t2,t1),t2_dom)

fwrite(promoter_anno[which(promoter_anno$V4 %in% t1_and_t2),c(1,2,3)] ,file = "8Cell_t1t2_promoter.bed", quote = F,sep="\t",col.names = F)
fwrite(promoter_anno[which(promoter_anno$V4 %in% t1_only),c(1,2,3)] ,file = "8Cell_t1_promoter.bed", quote = F,sep="\t",col.names = F)
fwrite(promoter_anno[which(promoter_anno$V4 %in% t2_only),c(1,2,3)] ,file = "8Cell_t2_promoter.bed", quote = F,sep="\t",col.names = F)

# Compare gene expr
dql <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Deng_Science_allstage_collapse.txt"),check.names = F)
tt_gene <- data.frame(name=c(t1_only,t2_only,t1_and_t2), class = c(rep("t1",length(t1_only)),rep("t2",length(t2_only)),rep("t1t2",length(t1_and_t2))))
tt_gene <- merge(tt_gene, dql, by.x = "name", by.y = "gene")[,c("name","class","8cell")]
ggboxplot(tt_gene,"class","8cell")


######morula
tf_promoter <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/standardize/3.tf_promoter/tf_promoter.tab"),check.names = F)
tf_promoter[,4:ncol(tf_promoter)] <- log2(tf_promoter[,4:ncol(tf_promoter)] + 1)
promoter_anno <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/1.dynamic_change/mm10_2k_0.5k_promoter.bed"))
tf_promoter <- merge(tf_promoter,promoter_anno, by.x=c("chr","start","end"), by.y=c("V1","V2","V3"))
t1_dom <- tf_promoter[which(tf_promoter$`Morula_Tfap2c` - tf_promoter$`Morula_Tead4` > 1),"V4"]
t2_dom <- tf_promoter[which(tf_promoter$`Morula_Tead4` - tf_promoter$`Morula_Tfap2c` > 1),"V4"]

promoter_anno <-  data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/fig1/1.dynamic_change/mm10_2k_0.5k_promoter.bed"))

t1 <- stage_tf_class[["Morula"]][["Tfap2c"]][["promoter"]]
t2 <- stage_tf_class[["Morula"]][["Tead4"]][["promoter"]]
t1_and_t2 <- intersect(t1,t2)
t1_only <- intersect(setdiff(t1,t2),t1_dom)
t2_only <- intersect(setdiff(t2,t1),t2_dom)

fwrite(promoter_anno[which(promoter_anno$V4 %in% t1_and_t2),c(1,2,3)] ,file = "Morula_t1t2_promoter.bed", quote = F,sep="\t",col.names = F)
fwrite(promoter_anno[which(promoter_anno$V4 %in% t1_only),c(1,2,3)] ,file = "Morula_t1_promoter.bed", quote = F,sep="\t",col.names = F)
fwrite(promoter_anno[which(promoter_anno$V4 %in% t2_only),c(1,2,3)] ,file = "Morula_t2_promoter.bed", quote = F,sep="\t",col.names = F)

# Compare gene expr
#dql <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Deng_Science_allstage_collapse.txt"),check.names = F)
tt_gene <- data.frame(name=c(t1_only,t2_only,t1_and_t2), class = c(rep("t1",length(t1_only)),rep("t2",length(t2_only)),rep("t1t2",length(t1_and_t2))))
tt_gene <- merge(tt_gene, dql, by.x = "name", by.y = "gene")[,c("name","class","earlyblast")]
ggboxplot(tt_gene,"class","earlyblast")


