#!/bin/sh
set -e
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
featureCounts -T 40  -s 0 -t exon -g gene_id -F GTF -a /media/wangjiachen/disk1/genome/sc/from_ens_Saccharomyces_cerevisiae.R64-1-1.57.gtf -o ${id}.count ${id}.bam
echo "count done"
sed '1d' ${id}.count | awk '{print $1"\t"$7}' > ${id}.rawcount
echo "rawcount done"
done
mkdir count ; mkdir rawcount ; mv *.count count ; mv *.rawcount rawcount ; mv *.summary count ; mkdir featurecount 
mv count featurecount ; mv rawcount featurecount ; cd featurecount ; cd count ; mkdir nouse_summary ; mv *.summary nouse_summary ; cd ../.. 
echo "featurecount done"
