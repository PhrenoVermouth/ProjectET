#!/bin/sh
THREAD=3
tempfifo=$$.fifo
trap 'exec 1000>&-;exec 1000<&+;rm -rf $tempfifo;exit 0' 2
mkfifo $tempfifo
exec 1000<>$tempfifo
for ((i=0;i<$THREAD;i++))
do
echo >&1000
done

wd=~/work_space/4.ProjectET/batch5_cr/3.align
cd $wd
file=`ls *sam`
for i in $file
do
read -u 1000
{
samtools sort -@ 10 -o ${i%.*}.sorted.bam $i
echo >&1000
}&
done
wait
echo "All Done!"
rm -rf $tempfifo

