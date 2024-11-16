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
/home/wangjiachen/anaconda3/envs/bowtie2/bin/bowtie2 -p 40 --no-mixed --no-discordant -x /media/wangjiachen/disk5/XWH/genome/Ensembl/Brachypodium_distachyon/bowtie2_index/bowtie2_index -1 ${id}_R1.fq.gz -2 ${id}_R2.fq.gz -S ${id}.sam
echo "map done"
cat ${id}.sam | grep -v 'Bd1_centromere_containing_Bradi1g41430' | grep -v 'KZ622971'| grep -v 'KZ622972' | grep -v 'KZ622973' | grep -v 'KZ622974' > ${id}_rmCMKZ.sam
echo "rmCMKZ done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -@ 40 -q 20  -bS ${id}_rmCMKZ.sam > ${id}_rmCMKZ_q20.bam
echo "q20 bam done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools sort -@ 40 -n ${id}_rmCMKZ_q20.bam -o ${id}_rmCMKZ_q20_s1.bam
echo "nsort1 done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools fixmate -@ 40 -m ${id}_rmCMKZ_q20_s1.bam ${id}_rmCMKZ_q20_s1_fixmate.bam
echo "fixmate done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools sort -@ 40 ${id}_rmCMKZ_q20_s1_fixmate.bam -o ${id}_rmCMKZ_q20_s1_fixmate_s2.bam
echo "sort2 done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools markdup -@ 40 -r ${id}_rmCMKZ_q20_s1_fixmate_s2.bam ${id}_rmCMKZ_q20_s1_fixmate_s2_markdup.bam
echo "markdup done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools index -@ 40 ${id}_rmCMKZ_q20_s1_fixmate_s2_markdup.bam
echo "index done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -@ 40 -F 4  -c ${id}.sam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -@ 40 -F 4  -c ${id}_rmCMKZ.sam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -@ 40 -F 4  -c ${id}_rmCMKZ_q20.bam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -@ 40 -F 4  -c ${id}_rmCMKZ_q20_s1_fixmate_s2_markdup.bam
echo "sam_count"
echo "rmCMKZ_count"
echo "q20_count"
echo "markdup_count"
rm ${id}_rmCMKZ_q20_s1_fixmate_s2.bam ; rm ${id}_rmCMKZ_q20_s1_fixmate.bam ; rm ${id}_rmCMKZ_q20_s1.bam ; rm ${id}_rmCMKZ_q20.bam ; rm ${id}_rmCMKZ.sam ; rm ${id}.sam
echo "findme"
done
mkdir Fastq4map ; mv *.fq.gz Fastq4map
echo "finalbam done"
#3.Normalization & BAM2BED
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
bamCoverage -e -p 40 -b ${id}.bam -o ${id}_RPKM_e.bw --binSize 10 --normalizeUsing RPKM
echo "RPKM bigwig done"
done
mkdir RPKM 
mv *.bw RPKM
echo "bigwig done"
echo  "All done"
