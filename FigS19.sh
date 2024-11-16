#Figure S19

#A
computeMatrix scale-regions -S \
WT_H3.bigwig \
ehd3_H3.bigwig \
sdg724_H3.bigwig \
ehd3sdg724_H3.bigwig \
-R Genes.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -p 30 -out H3.mat.gz
plotHeatmap -m H3.mat.gz -out H3_heatmap.pdf --sortRegions descend --sortUsingSamples 1 --colorList "white,#FF0000" --zMax 50


#B
plotProfile -m H3.mat.gz   -out H3_pattern.pdf  --perGroup --colors black red --samplesLabel NIP --legendLocation upper-right --yMin 0 --yMax 50 --plotHeight 11 --plotWidth 11

