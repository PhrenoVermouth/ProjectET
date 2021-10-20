#!/bin/sh
THREAD=4
tempfifo=$$.fifo
trap 'exec 1000>&-;exec 1000<&+;rm -rf $tempfifo;exit 0' 2
mkfifo $tempfifo
exec 1000<>$tempfifo
for ((i=0;i<$THREAD;i++))
do
echo >&1000
done

wd=~/work_space/4.ProjectET/batch9_cr/2.cutdata
nwd=~/work_space/4.ProjectET/batch9_cr/a.snp_align
mkdir -p $nwd
cd $wd
file=`ls *R1_cut.fq`
for i in $file
do
t=${i/R1/R2}
echo $i
echo $t
read -u 1000
{
bowtie2 -p 8 -x /home1/share/snpsplit/PWK_Phj_single_strain/PWK_Phj_bowtie2_index/PWK_Phj_N_masked  -t --no-mixed --no-discordant --no-unal -1 $i -2 $t -S $nwd/${i%%.R1*}.sam > $nwd/${i%%.R1*}.log2>&1
echo >&1000
}&
done
wait
echo "All Done!"
rm -rf $tempfifo
