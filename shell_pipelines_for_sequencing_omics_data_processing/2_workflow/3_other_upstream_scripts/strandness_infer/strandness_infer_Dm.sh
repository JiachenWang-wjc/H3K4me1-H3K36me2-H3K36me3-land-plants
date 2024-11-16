#!/bin/sh
set -e
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
infer_experiment.py -r /media/wangjiachen/disk6/test/scripts/genome/9_Dm/D_melanogaster_gene_col6.bed -i ${id}.bam
done
echo "strandness_infer done"

