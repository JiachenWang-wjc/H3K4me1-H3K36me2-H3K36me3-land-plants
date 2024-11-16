#!/bin/sh
set -e
#1.Quality control & trimming
mkdir QC
fastqc *.gz -o QC -t 40
ls -1 *.fastq.gz |awk -F ".fastq.gz" '{print $1}' | while read id
do
echo "${id}"
trim_galore --phred33 --fastqc -q 20  --stringency 3 --length 20 -o trimGalore_trim_s3_l20  ${id}.fastq.gz  -j 4
done
echo "trim done"
mkdir rawfastq ; mv *.gz rawfastq ; mv QC rawfastq ; cd trimGalore_trim_s3_l20 ; mv *.gz .. ; cd ..
ls -1 *_trimmed.fq.gz |awk -F "_trimmed.fq.gz" '{print $1}' | while read id
do
echo "${id}"
mv "${id}"_trimmed.fq.gz   "${id}".fq.gz
done
echo "rename done"
#2.Mapping & data cleaning
ls -1 *.fq.gz |awk -F ".fq.gz" '{print $1}' | while read id
do
echo "${id}"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/bowtie2 -p 40 -x /media/wangjiachen/disk5/XWH/genome/Ensembl/Brachypodium_distachyon/bowtie2_index/bowtie2_index -U ${id}.fq.gz -S ${id}.sam
echo "map done"
cat ${id}.sam | grep -v 'MtDNA' > ${id}_rmMtDNA.sam
echo "rmMtDNA done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -@ 40 -q 20 -bS ${id}_rmMtDNA.sam > ${id}_rmMtDNA_q20.bam
echo "q20 bam done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools sort -@ 40 ${id}_rmMtDNA_q20.bam -o ${id}_rmMtDNA_q20_s.bam
echo "sort done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools rmdup  -s ${id}_rmMtDNA_q20_s.bam ${id}_rmMtDNA_q20_s_rmdup.bam
echo "rmdup done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools index -@ 40 ${id}_rmMtDNA_q20_s_rmdup.bam
echo "index done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}.sam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}_rmMtDNA.sam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}_rmMtDNA_q20.bam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}_rmMtDNA_q20_s_rmdup.bam
echo "sam_count"
echo "rmMtDNA_count"
echo "q20_count"
echo "rmdup_count"
echo "reads_count_done"
rm ${id}_rmMtDNA_q20_s.bam ; rm ${id}_rmMtDNA_q20.bam ; rm ${id}_rmMtDNA.sam ; rm ${id}.sam
done
mkdir Fastq4map ; mv *.fq.gz Fastq4map
echo "finalbam done"
#3.Normalization & BAM2BED
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
	echo "${id}"
	bamCoverage -p 40 -b ${id}.bam -o ${id}_RPKM.bw --binSize 10 --normalizeUsing RPKM
	echo "RPKM bw done"
done
mkdir RPKM
mv *.bw RPKM
echo "bw done"
echo "All done"
