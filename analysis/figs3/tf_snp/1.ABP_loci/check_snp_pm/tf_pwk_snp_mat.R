########################################
# classify the TF peak into two classes, promoter and distal
# date: 2021.12.22 -- 12.23
# author: Jing Xiao
########################################


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
work_dir <- "D:/JLab/project/20211217_mouse_tfap2c_tead4/check_snp_pm"
setwd(work_dir)


# load packages ---------------------------------------------------------------
library(tidyverse)
library(stringr)


# load tf_PWK_detail.txt ------------------------------------------------------
tf_pwk_snp_detail_file <- list.files(pattern = "*PWK_detail.txt")
for (i in tf_pwk_snp_detail_file) {
    sample_name <- vapply(
        strsplit(i, "_"),
        function(x) paste(x[seq.int(2)], collapse = "_"),
        character(1L)
    )   # 去除第二个'_'后的内容
    detail_dta <- read_tsv(i, col_names = FALSE)
    colnames(detail_dta) <- c(
        "motif_chr", "motif_start", "motif_end",
        "snp_chr", "snp_start", "snp_end",
        "snp_index", "motif_region"
    )
    tf_pwk_snp_mat <- data.frame(
        matrix(
        0, 
        nrow = length(unique(detail_dta$motif_region)),
        ncol = 8
        )
    )
    colnames(tf_pwk_snp_mat) <- c(paste("index", 1:8, sep = "_"))
    rownames(tf_pwk_snp_mat) <- unique(detail_dta$motif_region)
    for (j in seq_len(dim(tf_pwk_snp_mat)[1])) {
        mat_row_name <- rownames(tf_pwk_snp_mat)[j]
        input_index <- which(detail_dta$motif_region == mat_row_name)
        snp_index <- pull(detail_dta[input_index, 7])
        tf_pwk_snp_mat[j, snp_index] <- 1
    }
    tf_pwk_snp_tib <-  tf_pwk_snp_mat %>%
        as_tibble() %>%
        mutate(motif_region = rownames(tf_pwk_snp_mat)) %>%
        relocate(motif_region, .before = "index_1")  

    write_tsv(
        tf_pwk_snp_tib,
        paste(sample_name, "SNP_tibble.txt", sep = "_"),
        col_names = TRUE
    )
}

###### count the SNP number in each position
# load tf_PWK_detail.txt ------------------------------------------------------
tf_pwk_snp_tibble_file <- list.files(pattern = "*SNP_tibble.txt")
position_stat_final <- c()
for (i in tf_pwk_snp_tibble_file) {
    allel_tf <- vapply(
        strsplit(i, "_"),
        function(x) paste(x[seq.int(1)], collapse = "_"),
        character(1L)
    )   # 去除第一个'_'后的内容
    allel <- str_sub(allel_tf, start = 1, end = 1)
    tf <- str_sub(allel_tf, start = 2)
    tibble_dta <- read_tsv(i, col_names = TRUE)
    tibble_df <- tibble_dta %>%
        select(starts_with("index")) %>%
        as.data.frame()
    column_sum <- colSums(tibble_df)
    position_stat_tibble <- tibble(
        tf = tf,
        allel = allel,
        position = seq_len(8),
        snp_num = column_sum
    ) %>%
    mutate(allelic_type = if_else(allel == "m", "maternal", "paternal")) %>%
    relocate(allelic_type, .after = "tf") %>%
    select(-allel)
    position_stat_final <- bind_rows(position_stat_final, position_stat_tibble)
}

position_stat_final <- position_stat_final %>%
    mutate(
        total_motif_num = rep(c(4322, 3382, 2314, 2813), each = 8),
        snp_percent = format(round(snp_num / total_motif_num, 3), nsmall = 3)
    ) %>%
    arrange(tf)
   
position_stat_final$allelic_type <- factor(
    position_stat_final$allelic_type,
    levels = c("paternal", "maternal")
)
position_stat_final$position <- as.character(position_stat_final$position)
position_stat_final$snp_percent <- as.numeric(position_stat_final$snp_percent)

write_tsv(position_stat_final, "TF_SNP_position_stat.txt", col_names = TRUE)

tf_motif_snp_ratio_bar <- ggplot(position_stat_final) +
	geom_bar(
		aes(x = position, y = snp_percent, fill = allelic_type),
		stat = "identity",
        position = "dodge"
	) +
	facet_wrap(. ~ tf)
pdf("tf_motif_snp_ratio_bar.pdf", width = 8, height = 8)
print(tf_motif_snp_ratio_bar)
dev.off()


###### count the SNP density in each motif, and then do the t.test()
# load tf_PWK_detail.txt ------------------------------------------------------
tf_pwk_snp_tibble_file <- list.files(pattern = "*SNP_tibble.txt")
motif_region_num <- c(4322, 3382, 2314, 2813)
snp_density_list <- list()
for (i in seq_len(4)) {
    allel_tf <- vapply(
        strsplit(tf_pwk_snp_tibble_file[i], "_"),
        function(x) paste(x[seq.int(1)], collapse = "_"),
        character(1L)
    )   # 去除第一个'_'后的内容
    tibble_dta <- read_tsv(tf_pwk_snp_tibble_file[i], col_names = TRUE)
    tibble_df <- tibble_dta %>%
        select(starts_with("index")) %>%
        as.data.frame()
    snp_density <- rowSums(tibble_df) / 8
    snp_density_new <- c(snp_density, rep(0, times = motif_region_num[i] - length(snp_density)))
    snp_density_list[[allel_tf]] <- snp_density_new
}

names(snp_density_list)
# [1] "mTead4"  "mTfap2c" "pTead4"  "pTfap2c"

tead4_wilcox_test <- wilcox.test(
    x = snp_density_list[["pTead4"]],
    y = snp_density_list[["mTead4"]],
    paired = FALSE
)
# 	Wilcoxon rank sum test with continuity correction

# data:  snp_density_list[["pTead4"]] and snp_density_list[["mTead4"]]
# W = 4734922, p-value = 2.262e-12
# alternative hypothesis: true location shift is not equal to 0

tfap2c_wilcox_test <- wilcox.test(
    x = snp_density_list[["pTfap2c"]],
    y = snp_density_list[["mTfap2c"]],
    paired = FALSE
)
# 	Wilcoxon rank sum test with continuity correction

# data:  snp_density_list[["pTfap2c"]] and snp_density_list[["mTfap2c"]]
# W = 4342835, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
