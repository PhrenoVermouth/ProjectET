########################################
# stat the number and percent of two TF allelic motif, group by with or without SNP
# date: 2021.12.22
# author: Jing Xiao
########################################


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
work_dir <- "D:/JLab/project/20211217_mouse_tfap2c_tead4/check_snp_pm"
setwd(work_dir)


# load packages ---------------------------------------------------------------
library(tidyverse)


# create stat data ------------------------------------------------------------
tf_motif_snp_stat_tib <- tibble(
	tf = c("Tead4", "Tead4", "Tfap2c", "Tfap2c"),
	allelic_type = c("paternal", "maternal", "paternal", "maternal"),
	total = c(2314, 4322, 2813, 3382),
	with_snp = c(141, 492, 161, 486)
)
tf_motif_snp_stat_tib <- tf_motif_snp_stat_tib %>%
	mutate(without_snp = total - with_snp) %>%
	gather(with_snp, without_snp, key = "snp_type", value = "num") %>%
	mutate(
		percent = num / total,
		anno_label = paste0(num, "/", total, ": ", format(round(percent, 3), nsmall = 3))
	) %>%
	arrange(tf, desc(allelic_type))
tf_motif_snp_stat_tib$allelic_type <- factor(
	tf_motif_snp_stat_tib$allelic_type,
	levels = c("paternal", "maternal")
)
write_tsv(
	tf_motif_snp_stat_tib, "tf_motif_snp_stat_tib.txt",
	col_names = TRUE, quote_escape = "none"
)


# anno_text <- data.frame(
#   x = c(1, 1, 2, 2, 1, 1, 2, 2),
#   y = c(0.98, 0.6, 0.98, 0.6, 0.98, 0.6, 0.98, 0.6),
#   label = tf_motif_snp_stat_tib$anno_label
# )
tf_motif_snp_stat_bar <- ggplot(tf_motif_snp_stat_tib) +
	geom_bar(
		aes(x = allelic_type, y = percent, fill = snp_type),
		stat = "identity"
	) +
	facet_wrap(. ~ tf)
pdf("tf_motif_snp_stat_bar.pdf", width = 8, height = 8)
print(tf_motif_snp_stat_bar)
dev.off()
