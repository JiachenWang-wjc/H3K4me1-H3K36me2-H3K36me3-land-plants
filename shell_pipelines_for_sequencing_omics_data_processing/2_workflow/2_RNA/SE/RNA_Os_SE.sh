#!/bin/sh
set -e
#1.quality control & trimming
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
#2.mapping & data cleaning
ls -1 *.fq.gz |awk -F ".fq.gz" '{print $1}' | while read id
do
	echo "${id}"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/hisat2 -p 40 -x /media/wangjiachen/disk1/genome/rice/MSU7_wjc/hisat2_index_withChrCM/hisat2_index_with_chrcm -U ${id}.fq.gz -S ${id}.sam
	echo "map done"
	cat ${id}.sam | grep -v 'ChrUn' | grep -v 'ChrSy'| grep -v 'ChrM' | grep -v 'ChrC' > ${id}_rmCMUnSy.sam
    echo "rmCMUnSy done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -@ 40 -q 20 -bS ${id}_rmCMUnSy.sam > ${id}_rmCMUnSy_q20.bam
	echo "q20 bam done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools sort -@ 40 ${id}_rmCMUnSy_q20.bam -o ${id}_rmCMUnSy_q20_s.bam
	echo "sort done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools rmdup -s ${id}_rmCMUnSy_q20_s.bam ${id}_rmCMUnSy_q20_s_rmdup.bam
	echo "rmdup done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools index -@ 40 ${id}_rmCMUnSy_q20_s_rmdup.bam
	echo "index done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}.sam
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmCMUnSy.sam
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmCMUnSy_q20.bam
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmCMUnSy_q20_s_rmdup.bam
	echo "sam_count"
	echo "rmCMUnSy_count"
	echo "q20_count"
	echo "rmdup_count"
	echo "reads_count_done"
	rm ${id}_rmCMUnSy_q20_s.bam ; rm ${id}_rmCMUnSy_q20.bam ; rm ${id}_rmCMUnSy.sam ; rm ${id}.sam
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



