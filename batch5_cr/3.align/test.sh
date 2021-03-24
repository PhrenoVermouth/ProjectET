#!/bin/sh
THREAD=5
tempfifo=$$.fifo
trap 'exec 1000>&-;exec 1000<&+;rm -rf $tempfifo;exit 0' 2
mkfifo $tempfifo
exec 1000<>$tempfifo
for ((i=0;i<$THREAD;i++))
do
echo >&1000
done

wd=~/work_space/4.ProjectET/batch5_cr/2.cutdata
nwd=~/work_space/4.ProjectET/batch5_cr/3.align
cd $wd
file=`ls *R1_cut.fq`
for i in $file
do
t=${i/R1/R2}
echo $i
echo $t
read -u 1000
{
nohup bowtie2 -p 12 -x /home1/share/bowtie2_index/mm10  --local --very-sensitive-local -t --no-mixed --no-discordant --no-unal -1 $i -2 $t -S $nwd/${i%%.*}.sam > $nwd/${i%%.*}.log 
echo >&1000
}&
done
wait
echo "All Done!"
rm -rf $tempfifo

