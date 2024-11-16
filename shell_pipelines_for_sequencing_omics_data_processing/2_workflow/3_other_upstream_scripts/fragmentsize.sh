#!/bin/sh
set -e
ls -1 *.bam > bam_name.txt
ls -1 *.bam | awk -F "_rm" '{print $1}' >sample_name.txt
sed -i ':label;N;s/\n/ /;b label'  bam_name.txt
sed -i ':label;N;s/\n/ /;b label'  sample_name.txt
cat bam_name.txt |  while read id
do
echo "${id}"
cat sample_name.txt |  while read sample
do
echo "${sample}"
bamPEFragmentSize  --bamfiles ${id} --histogram all_bam_fragsize.pdf  --numberOfProcessors 4 --maxFragmentLength 500 --samplesLabel ${sample}
done
done
mkdir fragmentsize ; mv *.pdf fragmentsize ; mv bam_name.txt sample_name.txt fragmentsize
echo "fragmentsize done" 