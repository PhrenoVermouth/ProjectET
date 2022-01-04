########## 2021.12.22
# check the frequency of PWK SNP in Tfap2c / Tead4 paternal or maternal biased motif region
# bedtools intersect
# bedtools intersect, version: v2.28.0

work_dir="/home1/gyang/work_space/4.ProjectET/analysis/figs3/tf_snp/1.ABP_loci"
snp_dir="/home1/share/snpsplit/PWK_Phj_single_strain/all_SNPs_PWK_PhJ_GRCm38.txt.gz"
cd ${work_dir}

mkdir check_snp_pm

# get the PWK SNP region
nohup zcat /home1/share/snpsplit/PWK_Phj_single_strain/all_SNPs_PWK_PhJ_GRCm38.txt.gz | awk 'BEGIN { FS = OFS = "\t" } \
{ print $2, $3, $3 }' > ${work_dir}/check_snp_pm/PWK_SNP_mm10.bed &

# retain the motif region whose lenght is 8
ls *_abs.bed | while read id;
do
    sample=${id%%.*}
    nohup awk 'BEGIN { FS = OFS = "\t" } { $4 = $3 - $2 + 1 } { if ($4 == 8) print $1, $2, $3 }' \
    ${id} > ${work_dir}/check_snp_pm/${sample}_filter.bed &
done

cd ${work_dir}/check_snp_pm
ls *_filter.bed | while read id;
do
    sample=${id%%_*}
    nohup bedtools intersect -a ${id} -b PWK_SNP_mm10.bed -wa -wb > ${sample}_PWK_info.txt &
    nohup bedtools intersect -a ${id} -b PWK_SNP_mm10.bed -wa -u > ${sample}_PWK_uni.txt &
done

# motif region的数量
wc -l *_filter.bed
#   4322 mTead4_abs_filter.bed
#   3382 mTfap2c_abs_filter.bed
#   2314 pTead4_abs_filter.bed
#   2813 pTfap2c_abs_filter.bed

# 含有SNP的motif region的数量
wc -l *_PWK_uni.txt
#   492 mTead4_PWK_uni.txt
#   486 mTfap2c_PWK_uni.txt
#   141 pTead4_PWK_uni.txt
#   161 pTfap2c_PWK_uni.txt

# 获取SNP在motif region的index（即位于第几个位置）
ls *_PWK_info.txt | while read id;
do
    sample=${id%_*}
    nohup awk 'BEGIN { FS = OFS = "\t" } { $7 = $5 - $2 + 1; $8 = $1"_"$2"_"$3; print $0 }' ${id} \
    > ${sample}_detail.txt &
done

# create tf_pwk_snp_mat.R
touch tf_pwk_snp_mat.R
nohup Rscript4 tf_pwk_snp_mat.R &
