#!/bin/sh  
set -e  

mkdir -p bam2bed  

for bamfile in *_rm*.bam; do  
    id=${bamfile%_rm*.bam} 
    echo "${id}"  
    /home/wangjiachen/anaconda3/envs/bedtools/bin/bedtools bamtobed -i "${bamfile}" > "${id}.bed"  
    echo "bam2bed done"    
    mv "${id}.bed" bam2bed/  
done  

echo "done"  

