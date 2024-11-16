#!/bin/sh
set -e
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
infer_experiment.py -r /media/wangjiachen/disk1/genome/rice/MSU7_wjc/rice_gene_noUnSy.bed -i ${id}.bam
done
echo "strandness_infer done"
