#!/bin/sh
set -e
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
infer_experiment.py -r /media/wangjiachen/disk1/genome/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed -i ${id}.bam
done
echo "strandness_infer done"
