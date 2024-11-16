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
	/home/wangjiachen/anaconda3/envs/hisat2/bin/hisat2 -p 40 -x /media/wangjiachen/disk1/genome/Cr/hisat2_Cr_lz_index/hisat2_index_Cr -U ${id}.fq.gz -S ${id}.sam
	echo "map done"
	cat ${id}.sam | grep -v 'KZ454945' | grep -v 'KZ454946' | grep -v 'KZ454947' | grep -v 'KZ454948' | grep -v 'KZ454949' | grep -v 'KZ454950' | grep -v 'KZ454951' | grep -v 'KZ454952' | grep -v 'KZ454953' | grep -v 'KZ454954' | grep -v 'KZ454955' | grep -v 'KZ454956' | grep -v 'KZ454957' | grep -v 'KZ454958' | grep -v 'KZ454959' | grep -v 'KZ454960' | grep -v 'KZ454961' | grep -v 'KZ454962' | grep -v 'KZ454963' | grep -v 'KZ454964' | grep -v 'KZ454965' | grep -v 'KZ454966' | grep -v 'KZ454967' | grep -v 'KZ454968' | grep -v 'KZ454969' | grep -v 'KZ454970' | grep -v 'KZ454971' | grep -v 'KZ454972' | grep -v 'KZ454973' | grep -v 'KZ454974' | grep -v 'KZ454975' | grep -v 'KZ454976' | grep -v 'KZ454977' | grep -v 'KZ454978' | grep -v 'KZ454979' | grep -v 'KZ454980' > ${id}_rmKZ4549XX.sam
    echo "rmKZ4549XX done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -@ 40 -q 20 -bS ${id}_rmKZ4549XX.sam > ${id}_rmKZ4549XX_q20.bam
	echo "q20 bam done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools sort -@ 40 ${id}_rmKZ4549XX_q20.bam -o ${id}_rmKZ4549XX_q20_s.bam
	echo "sort done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools rmdup -s ${id}_rmKZ4549XX_q20_s.bam ${id}_rmKZ4549XX_q20_s_rmdup.bam
	echo "rmdup done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools index -@ 40 ${id}_rmKZ4549XX_q20_s_rmdup.bam
	echo "index done"
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}.sam
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmKZ4549XX.sam
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmKZ4549XX_q20.bam
	/home/wangjiachen/anaconda3/envs/hisat2/bin/samtools view -F 4 -@ 40 -c ${id}_rmKZ4549XX_q20_s_rmdup.bam
	echo "sam_count"
	echo "rmKZ4549XX_count"
	echo "q20_count"
	echo "rmdup_count"
	echo "reads_count_done"
	rm ${id}_rmKZ4549XX_q20_s.bam ; rm ${id}_rmKZ4549XX_q20.bam ; rm ${id}_rmKZ4549XX.sam ; rm ${id}.sam
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

