#!/bin/sh
set -e
#1.Quality control & trimming
mkdir QC
fastqc *.gz -o QC -t 40
ls -1 *_R1.fastq.gz |awk -F "_R1.fastq.gz" '{print $1}' | while read id
do
echo "${id}"
trim_galore --phred33 --fastqc -q 20  --stringency 3 --length 20 -o trimGalore_trim_s3_l20 --paired ${id}_R1.fastq.gz ${id}_R2.fastq.gz -j 4
done
echo "trim done"
mkdir rawfastq ; mv *.gz rawfastq ; mv QC rawfastq ; cd trimGalore_trim_s3_l20 ; mv *.gz .. ; cd ..
ls -1 *_R1_val_1.fq.gz |awk -F "_R1_val_1.fq.gz" '{print $1}' | while read id
do
echo "${id}"
mv "${id}"_R1_val_1.fq.gz "${id}"_R1.fq.gz
mv "${id}"_R2_val_2.fq.gz "${id}"_R2.fq.gz
done
echo "rename done"
#2.Mapping & data cleaning
ls -1 *_R1.fq.gz |awk -F "_R1.fq.gz" '{print $1}' | while read id
do
echo "${id}"
/home/wangjiachen/anaconda3/envs/hisat2/bin/hisat2 -p 40 --dta --rna-strandness RF -x /media/wangjiachen/disk5/XWH/genome/Ensembl/D_melanogaster/hisat2_index/hisat2_index -1 ${id}_R1.fq.gz -2 ${id}_R2.fq.gz -S ${id}.sam
echo "map done"
cat ${id}.sam | grep -v 'mitochondrion_genome' > ${id}_rmMt.sam
echo "rmMt done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -@ 40 -q 20 -bS ${id}_rmMt.sam > ${id}_rmMt_q20.bam
echo "q20 bam done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools sort -@ 40 -n ${id}_rmMt_q20.bam -o ${id}_rmMt_q20_s1.bam
echo "nsort1 done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools fixmate -@ 40 -m ${id}_rmMt_q20_s1.bam ${id}_rmMt_q20_s1_fixmate.bam
echo "fixmate done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools sort -@ 40 ${id}_rmMt_q20_s1_fixmate.bam -o ${id}_rmMt_q20_s1_fixmate_s2.bam
echo "sort2 done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools markdup -@ 40 -r ${id}_rmMt_q20_s1_fixmate_s2.bam ${id}_rmMt_q20_s1_fixmate_s2_markdup.bam
echo "markdup done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools sort -@ 40 ${id}_rmMt_q20.bam -o ${id}_rmMt_q20_s.bam
echo "sort done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools index -@ 40 ${id}_rmMt_q20_s.bam
echo "index done"
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}.sam
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmMt.sam
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmMt_q20.bam
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmMt_q20_s.bam
/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmMt_q20_s1_fixmate_s2_markdup.bam
echo "reads count done"
rm ${id}.sam ; rm ${id}_rmMt.sam ; rm ${id}_rmMt_q20.bam ; rm ${id}_rmMt_q20_s1_fixmate_s2_markdup.bam
rm ${id}_rmMt_q20_s1.bam ; rm ${id}_rmMt_q20_s1_fixmate.bam ; rm ${id}_rmMt_q20_s1_fixmate_s2.bam
echo "sam_count"
echo "rmMt_count"
echo "q20_count"
echo "sort_count"
echo "rmdup_count"
echo "findme"
done
mkdir Fastq4map ; mv *.fq.gz Fastq4map
echo "finalbam done"
#3.Normalization
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
bamCoverage -p 40 -b ${id}.bam -o ${id}_RPKM.bw --binSize 10 --normalizeUsing RPKM
echo "RPKM bigwig done"
done
mkdir RPKM ; mv *.bw RPKM
echo "bigwig done"
echo "All done"

