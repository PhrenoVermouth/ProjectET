setwd("~/work_space/4.ProjectET/analysis/figs5")

library(ggpubr)
library(ggrepel)
library(data.table)

###########  By loop
# 28037 enhancer.bed
# 10595 Tead4_enhancer.bed
# 5711 Tfap2c_enhancer.bed
# 4451 overlap
### By gene ($4+sort+uniq)
# 7743 enhancer.bed
# 5509 Tead4_enhancer.bed
# 3677 Tfap2c_enhancer.bed
# 3102 overlap
# df_e <- data.frame(class=c("None","Tfap2c","Tead4"), num=c(28037-10595-5711, 5711,10595))
# df_e$perc <- round(df_e$num/28037*100,2)
# labs <- paste0(df_e$class, " (", df_e$perc, "%)")
# ggpie(df_e, "num", label = labs, lab.pos = "in",
#       fill = "class", color = "white",
#       palette = c("grey70", "Salmon","RosyBrown1") )
# ggsave("TT_occupancy_enhancer_by_enhancer.pdf",dpi=300,width = 4,height = 6)
#
#
# pdf("Venn_TT_occupancy_enhancer_by_enhancer.pdf")
# D1<-venn.diagram(list(ESC=esc_pks,TSC=tsc_pks),filename=NULL,lwd=1,lty=1,col=c('red','blue'),
#                  fill=c('red','blue'),cat.col=c('red','blue'),rotation.degree=0,main = "Distal peaks shared by ESC and TSC")
# grid.draw(D1)
# dev.off()

###########  Loop repeats

setwd("~/work_space/4.ProjectET/analysis/figs7/1.Venn_loop/results_repeats_subfamily")


all_loop_rep <- list.files(pattern = "?subcounts$")

# stage_tf_class <- list(Eight_cell = list(Tfap2c=list(promoter=c(),distal=c()), Tead4=list(promoter=c(),distal=c())),
#                        Morula = list(Tfap2c=list(promoter=c(),distal=c()), Tead4=list(promoter=c(),distal=c())))

rep_counts = data.frame(fread("None_repeats.subcounts"))
colnames(rep_counts) <- c("None_num","class")


for (i in all_loop_rep[2:length(all_loop_rep)]){
  temp=data.frame(fread(i))
  colnames(temp) <- c(paste0(gsub("_repeats.subcounts","",i),"_num"),"class")
  rep_counts <- merge(rep_counts,temp, by="class",all = T)
}
rep_counts[is.na(rep_counts)] <- 0
#rep_counts$rsum <- rowSums(rep_counts[,2:9])
rep_counts$Nonefc <- log2(rep_counts$None_num +1) - log2(rep_counts$random_None_num +1)
rep_counts$Tead4fc <- log2(rep_counts$Tead4_specific_num +1) - log2(rep_counts$random_Tead4_specific_num +1)
rep_counts$Tfap2cfc <- log2(rep_counts$Tfap2c_specific_num +1) - log2(rep_counts$random_Tfap2c_specific_num +1)
rep_counts$TTfc <- log2(rep_counts$TT_intersect_num +1) - log2(rep_counts$random_TT_intersect_num +1)

# [1] "class"                      "None_num"                   "random_None_num"            "random_Tead4_specific_num"
# [5] "random_Tfap2c_specific_num" "random_TT_intersect_num"    "Tead4_specific_num"         "Tfap2c_specific_num"
# [9] "TT_intersect_num"           "Tfap2cfc"                   "TTfc"                       "Nonefc"
# [13] "Tead4fc"
options(ggrepel.max.overlaps = Inf)

plt <- rep_counts[,c("class", "None_num" ,"Nonefc" )]
plt$label <- ""
plt$label[which(plt$None_num>100 & plt$Nonefc > 2)] <- plt$class[which(plt$None_num>100 & plt$Nonefc > 2)]
plt$None_num <- log2(plt$None_num + 1)

p <-
  ggscatter(plt, x = "None_num", y = "Nonefc",alpha = 1,color = ifelse(plt$label == "", "grey50", "red"), palette = c("Firebrick4","grey"),label = "label",repel = T, legend="none" ) +
  geom_hline(linetype = "dashed",yintercept = 2,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = 6.66,color = "grey") +
  xlim(-0.5,12) + ylim(-7,7)+coord_cartesian(expand = F)
p + border()
ggsave("None_loop_enhancer.pdf", width = 6, height = 3.5)


# [1] "class"                      "None_num"                   "random_None_num"            "random_Tead4_specific_num"
# [5] "random_Tfap2c_specific_num" "random_TT_intersect_num"    "Tead4_specific_num"         "Tfap2c_specific_num"
# [9] "TT_intersect_num"           "Tfap2cfc"                   "TTfc"                       "Nonefc"
# [13] "Tead4fc"

plt <- rep_counts[,c("class", "Tfap2c_specific_num"  ,"Tfap2cfc" )]
plt$label <- ""
plt$label[which(plt$Tfap2c_specific_num>100 & plt$Tfap2cfc > 2)] <- plt$class[which(plt$Tfap2c_specific_num>100 & plt$Tfap2cfc > 2)]
plt$Tfap2c_specific_num <- log2(plt$Tfap2c_specific_num + 1)

p <-
  ggscatter(plt, x = "Tfap2c_specific_num", y = "Tfap2cfc",alpha = 1,color = ifelse(plt$label == "", "grey50", "red"), palette = c("Firebrick4","grey"),label = "label",repel = T, legend="none" ) +
  geom_hline(linetype = "dashed",yintercept = 2,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = 6.66,color = "grey") +
  xlim(-0.5,8) + ylim(-7,7)+coord_cartesian(expand = F)
p + border()
ggsave("Tfap2c_loop_enhancer.pdf", width = 6, height = 3.5)


plt <- rep_counts[,c("class", "Tead4_specific_num"  ,"Tead4fc" )]
plt$label <- ""
plt$label[which(plt$Tead4_specific_num>100 & plt$Tead4fc > 2)] <- plt$class[which(plt$Tead4_specific_num>100 & plt$Tead4fc > 2)]
plt$Tead4_specific_num <- log2(plt$Tead4_specific_num + 1)

p <-
  ggscatter(plt, x = "Tead4_specific_num", y = "Tead4fc",alpha = 1,color = ifelse(plt$label == "", "grey50", "red"), palette = c("Firebrick4","grey"),label = "label",repel = T, legend="none" ) +
  geom_hline(linetype = "dashed",yintercept = 2,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = 6.66,color = "grey") +
  xlim(-0.5,12) + ylim(-7,7)+coord_cartesian(expand = F)
p + border()
ggsave("Tead4_loop_enhancer.pdf", width = 6, height = 3.5)


plt <- rep_counts[,c("class", "TT_intersect_num"  ,"TTfc" )]
plt$label <- ""
plt$label[which(plt$TT_intersect_num>100 & plt$TTfc > 2)] <- plt$class[which(plt$TT_intersect_num>100 & plt$TTfc > 2)]
plt$TT_intersect_num <- log2(plt$TT_intersect_num + 1)

p <-
  ggscatter(plt, x = "TT_intersect_num", y = "TTfc",alpha = 1,color = ifelse(plt$label == "", "grey50", "red"), palette = c("Firebrick4","grey"),label = "label",repel = T, legend="none" ) +
  geom_hline(linetype = "dashed",yintercept = 2,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = 6.66,color = "grey") +
  xlim(-0.5,11) + ylim(-7,7)+coord_cartesian(expand = F)
p + border()
ggsave("TT_loop_enhancer.pdf", width = 6, height = 3.5)

#
# ggplot(plt, aes(None_num, Nonefc, label = label)) +
#   geom_text_repel() +
#   geom_point(color = ifelse(plt$label == "", "grey50", "red"))
