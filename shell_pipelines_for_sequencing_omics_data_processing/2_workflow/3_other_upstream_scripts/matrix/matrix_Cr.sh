#!/bin/sh
set -e
cd RPKM
ls *.bw >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bw | awk -F "_rm" '{print $1}' >bw_name2.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name2.txt
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name2.txt |  while read slabel
do
echo "${slabel}"
computeMatrix scale-regions \
-S ${id} \
-p 40 -R /media/wangjiachen/disk1/genome/Cr/Cr_gene.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -out all_matrix.mat.gz
echo "matrix1 done"
plotHeatmap -m all_matrix.mat.gz -out all_matrix_heatmap.pdf --sortRegions keep --colorList "white,#77AC30"
echo "plot1 done"
plotProfile -m all_matrix.mat.gz   -out all_matrix_pattern.pdf  \
--perGroup --legendLocation upper-right  \
--colors grey black blue pink darkred red yellow green purple brown olive orange "#FF00FF" "#CCFB5D" "#46C7C7" "#000080" "#56A5EC" \
--samplesLabel ${slabel}
done
done
echo "plot2 done"
cd ../
mkdir allmatrix ; mv RPKM/*.pdf allmatrix ; mv RPKM/all_matrix.mat.gz allmatrix
