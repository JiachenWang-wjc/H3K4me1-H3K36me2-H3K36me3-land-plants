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
/home/wangjiachen/anaconda3/envs/bowtie2/bin/bowtie2 -p 40 -x /media/wangjiachen/disk1/genome/Pp/bowtie2_Pp_tx/bowtie2_index -U ${id}.fq.gz -S ${id}.sam
echo "map done"
cat ${id}.sam | grep -v 'scaffold_37' | grep -v 'scaffold_38' | grep -v 'scaffold_50' | grep -v 'scaffold_51' | grep -v 'scaffold_47' | grep -v 'scaffold_48' | grep -v 'scaffold_49' | grep -v 'scaffold_56' | grep -v 'scaffold_52' | grep -v 'scaffold_53' | grep -v 'scaffold_62' | grep -v 'scaffold_64' | grep -v 'scaffold_67' | grep -v 'scaffold_76' | grep -v 'scaffold_173' | grep -v 'scaffold_165' | grep -v 'scaffold_79' | grep -v 'scaffold_81' | grep -v 'scaffold_85' | grep -v 'scaffold_105' | grep -v 'scaffold_89' | grep -v 'scaffold_93' | grep -v 'scaffold_95' | grep -v 'scaffold_97' | grep -v 'scaffold_876' | grep -v 'scaffold_100' | grep -v 'scaffold_101' | grep -v 'scaffold_104' | grep -v 'scaffold_111' | grep -v 'scaffold_113' | grep -v 'scaffold_115' | grep -v 'scaffold_116' | grep -v 'scaffold_117' | grep -v 'scaffold_121' | grep -v 'scaffold_123' | grep -v 'scaffold_124' | grep -v 'scaffold_127' | grep -v 'scaffold_128' | grep -v 'scaffold_135' | grep -v 'scaffold_137' | grep -v 'scaffold_142' | grep -v 'scaffold_145' | grep -v 'scaffold_149' | grep -v 'scaffold_152' | grep -v 'scaffold_162' | grep -v 'scaffold_176' | grep -v 'scaffold_177' | grep -v 'scaffold_179' | grep -v 'scaffold_180' | grep -v 'scaffold_187' | grep -v 'scaffold_188' | grep -v 'scaffold_194' | grep -v 'scaffold_200' | grep -v 'scaffold_202' | grep -v 'scaffold_205' | grep -v 'scaffold_203' | grep -v 'scaffold_214' | grep -v 'scaffold_215' | grep -v 'scaffold_441' | grep -v 'scaffold_222' | grep -v 'scaffold_224' | grep -v 'scaffold_226' | grep -v 'scaffold_230' | grep -v 'scaffold_241' | grep -v 'scaffold_228' | grep -v 'scaffold_229' | grep -v 'scaffold_234' | grep -v 'scaffold_240' | grep -v 'scaffold_242' | grep -v 'scaffold_246' | grep -v 'scaffold_244' | grep -v 'scaffold_245' | grep -v 'scaffold_261' | grep -v 'scaffold_254' | grep -v 'scaffold_255' | grep -v 'scaffold_256' | grep -v 'scaffold_267' | grep -v 'scaffold_259' | grep -v 'scaffold_270' | grep -v 'scaffold_271' | grep -v 'scaffold_272' | grep -v 'scaffold_274' | grep -v 'scaffold_1175' | grep -v 'scaffold_278' | grep -v 'scaffold_281' | grep -v 'scaffold_292' | grep -v 'scaffold_293' | grep -v 'scaffold_298' | grep -v 'scaffold_302' | grep -v 'scaffold_313' | grep -v 'scaffold_310' | grep -v 'scaffold_312' | grep -v 'scaffold_317' | grep -v 'scaffold_318' | grep -v 'scaffold_329' | grep -v 'scaffold_332' | grep -v 'scaffold_421' | grep -v 'scaffold_334' | grep -v 'scaffold_338' | grep -v 'scaffold_342' | grep -v 'scaffold_346' | grep -v 'scaffold_349' | grep -v 'scaffold_350' | grep -v 'scaffold_352' | grep -v 'scaffold_353' | grep -v 'scaffold_354' | grep -v 'scaffold_355' | grep -v 'scaffold_361' | grep -v 'scaffold_358' | grep -v 'scaffold_360' | grep -v 'scaffold_366' | grep -v 'scaffold_365' | grep -v 'scaffold_370' | grep -v 'scaffold_447' | grep -v 'scaffold_374' | grep -v 'scaffold_375' | grep -v 'scaffold_384' | grep -v 'scaffold_380' | grep -v 'scaffold_383' | grep -v 'scaffold_507' | grep -v 'scaffold_435' | grep -v 'scaffold_391' | grep -v 'scaffold_392' | grep -v 'scaffold_395' | grep -v 'scaffold_396' | grep -v 'scaffold_397' | grep -v 'scaffold_398' | grep -v 'scaffold_427' | grep -v 'scaffold_403' | grep -v 'scaffold_419' | grep -v 'scaffold_404' | grep -v 'scaffold_407' | grep -v 'scaffold_410' | grep -v 'scaffold_418' | grep -v 'scaffold_429' | grep -v 'scaffold_422' | grep -v 'scaffold_423' | grep -v 'scaffold_425' | grep -v 'scaffold_456' | grep -v 'scaffold_430' | grep -v 'scaffold_433' | grep -v 'scaffold_442' | grep -v 'scaffold_468' | grep -v 'scaffold_437' | grep -v 'scaffold_439' | grep -v 'scaffold_450' | grep -v 'scaffold_443' | grep -v 'scaffold_446' | grep -v 'scaffold_453' | grep -v 'scaffold_454' | grep -v 'scaffold_455' | grep -v 'scaffold_537' | grep -v 'scaffold_457' | grep -v 'scaffold_458' | grep -v 'scaffold_470' | grep -v 'scaffold_459' | grep -v 'scaffold_464' | grep -v 'scaffold_466' | grep -v 'scaffold_471' | grep -v 'scaffold_472' | grep -v 'scaffold_474' | grep -v 'scaffold_477' | grep -v 'scaffold_480' | grep -v 'scaffold_482' | grep -v 'scaffold_535' | grep -v 'scaffold_489' | grep -v 'scaffold_519' | grep -v 'scaffold_493' | grep -v 'scaffold_612' | grep -v 'scaffold_499' | grep -v 'scaffold_501' | grep -v 'scaffold_504' | grep -v 'scaffold_506' | grep -v 'scaffold_511' | grep -v 'scaffold_1024' | grep -v 'scaffold_518' | grep -v 'scaffold_573' | grep -v 'scaffold_523' | grep -v 'scaffold_524' | grep -v 'scaffold_545' | grep -v 'scaffold_527' | grep -v 'scaffold_530' | grep -v 'scaffold_531' | grep -v 'scaffold_533' | grep -v 'scaffold_542' | grep -v 'scaffold_546' | grep -v 'scaffold_547' | grep -v 'scaffold_548' | grep -v 'scaffold_569' | grep -v 'scaffold_550' | grep -v 'scaffold_592' | grep -v 'scaffold_554' | grep -v 'scaffold_564' | grep -v 'scaffold_565' | grep -v 'scaffold_567' | grep -v 'scaffold_568' | grep -v 'scaffold_571' | grep -v 'scaffold_844' | grep -v 'scaffold_710' | grep -v 'scaffold_576' | grep -v 'scaffold_578' | grep -v 'scaffold_581' | grep -v 'scaffold_755' | grep -v 'scaffold_585' | grep -v 'scaffold_589' | grep -v 'scaffold_614' | grep -v 'scaffold_593' | grep -v 'scaffold_866' | grep -v 'scaffold_868' | grep -v 'scaffold_598' | grep -v 'scaffold_613' | grep -v 'scaffold_616' | grep -v 'scaffold_625' | grep -v 'scaffold_635' | grep -v 'scaffold_637' | grep -v 'scaffold_638' | grep -v 'scaffold_693' | grep -v 'scaffold_651' | grep -v 'scaffold_858' | grep -v 'scaffold_656' | grep -v 'scaffold_660' | grep -v 'scaffold_669' | grep -v 'scaffold_674' | grep -v 'scaffold_675' | grep -v 'scaffold_677' | grep -v 'scaffold_810' | grep -v 'scaffold_689' | grep -v 'scaffold_691' | grep -v 'scaffold_724' | grep -v 'scaffold_695' | grep -v 'scaffold_699' | grep -v 'scaffold_712' | grep -v 'scaffold_714' | grep -v 'scaffold_738' | grep -v 'scaffold_1005' | grep -v 'scaffold_752' | grep -v 'scaffold_754' | grep -v 'scaffold_760' | grep -v 'scaffold_769' | grep -v 'scaffold_808' | grep -v 'scaffold_772' | grep -v 'scaffold_775' | grep -v 'scaffold_776' | grep -v 'scaffold_778' | grep -v 'scaffold_783' | grep -v 'scaffold_791' | grep -v 'scaffold_792' | grep -v 'scaffold_793' | grep -v 'scaffold_794' | grep -v 'scaffold_1160' | grep -v 'scaffold_804' | grep -v 'scaffold_816' | grep -v 'scaffold_818' | grep -v 'scaffold_829' | grep -v 'scaffold_831' | grep -v 'scaffold_841' | grep -v 'scaffold_907' | grep -v 'scaffold_851' | grep -v 'scaffold_912' | grep -v 'scaffold_909' | grep -v 'scaffold_877' | grep -v 'scaffold_920' | grep -v 'scaffold_881' | grep -v 'scaffold_995' | grep -v 'scaffold_893' | grep -v 'scaffold_898' | grep -v 'scaffold_1004' | grep -v 'scaffold_900' | grep -v 'scaffold_930' | grep -v 'scaffold_903' | grep -v 'scaffold_927' | grep -v 'scaffold_949' | grep -v 'scaffold_959' | grep -v 'scaffold_918' | grep -v 'scaffold_1015' | grep -v 'scaffold_942' | grep -v 'scaffold_1127' | grep -v 'scaffold_970' | grep -v 'scaffold_1140' | grep -v 'scaffold_937' | grep -v 'scaffold_1097' | grep -v 'scaffold_984' | grep -v 'scaffold_944' | grep -v 'scaffold_982' | grep -v 'scaffold_1033' | grep -v 'scaffold_960' | grep -v 'scaffold_1051' | grep -v 'scaffold_1119' | grep -v 'scaffold_1067' | grep -v 'scaffold_967' | grep -v 'scaffold_993' | grep -v 'scaffold_1053' | grep -v 'scaffold_1035' | grep -v 'scaffold_1060' | grep -v 'scaffold_1070' | grep -v 'scaffold_1088' | grep -v 'scaffold_1177' | grep -v 'scaffold_991' | grep -v 'scaffold_1105' | grep -v 'scaffold_1034' | grep -v 'scaffold_1150' | grep -v 'scaffold_1006' | grep -v 'scaffold_1121' | grep -v 'scaffold_1137' | grep -v 'scaffold_1016' | grep -v 'scaffold_1083' | grep -v 'scaffold_1078' | grep -v 'scaffold_1073' | grep -v 'scaffold_1045' | grep -v 'scaffold_1129' | grep -v 'scaffold_1090' | grep -v 'scaffold_1155' | grep -v 'scaffold_1151' | grep -v 'scaffold_1080' | grep -v 'scaffold_1167' | grep -v 'scaffold_1126' | grep -v 'scaffold_1172' | grep -v 'scaffold_1164' | grep -v 'scaffold_1176' | grep -v 'scaffold_1147' | grep -v 'scaffold_1153' | grep -v 'scaffold_1154' | grep -v 'scaffold_1157' | grep -v 'scaffold_1168' | grep -v 'scaffold_1173' | grep -v 'scaffold_1174' | grep -v 'scaffold_1178' | grep -v 'scaffold_1180' | grep -v 'scaffold_1181' | grep -v 'scaffold_1182'> ${id}_rmscaffold.sam
echo "rmscaffold done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -@ 40 -q 20 -bS ${id}_rmscaffold.sam > ${id}_rmscaffold_q20.bam
echo "q20 bam done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools sort -@ 40 ${id}_rmscaffold_q20.bam -o ${id}_rmscaffold_q20_s.bam
echo "sort done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools rmdup  -s ${id}_rmscaffold_q20_s.bam ${id}_rmscaffold_q20_s_rmdup.bam
echo "rmdup done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools index -@ 40 ${id}_rmscaffold_q20_s_rmdup.bam
echo "index done"
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}.sam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}_rmscaffold.sam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}_rmscaffold_q20.bam
/home/wangjiachen/anaconda3/envs/bowtie2/bin/samtools view -F 4 -@ 40 -c ${id}_rmscaffold_q20_s_rmdup.bam
echo "sam_count"
echo "rmscaffold_count"
echo "q20_count"
echo "rmdup_count"
echo "reads_count_done"
rm ${id}_rmscaffold_q20_s.bam ; rm ${id}_rmscaffold_q20.bam ; rm ${id}_rmscaffold.sam ; rm ${id}.sam
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
