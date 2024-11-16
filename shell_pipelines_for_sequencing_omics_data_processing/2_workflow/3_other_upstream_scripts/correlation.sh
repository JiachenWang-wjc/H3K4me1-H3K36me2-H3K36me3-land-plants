#!/bin/sh
set -e
/home/wangjiachen/anaconda3/envs/deeptools2/bin/multiBamSummary bins -bs 1000 -p 40 --bamfiles *.bam -out readCounts1.npz --outRawCounts readCounts1.tab
echo "corrlation matrix done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_PearsonCorr1.png --outFileCorMatrix PearsonCorr1.tab
echo "1 done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_PearsonCorr2.png   --outFileCorMatrix PearsonCorr2.tab
echo "2 done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_SpearmanCorr1.png --outFileCorMatrix SpearmanCorr1.tab
echo "3 done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_SpearmanCorr2.png --outFileCorMatrix SpearmanCorr2.tab
echo "4 done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_PearsonCorr3.png --outFileCorMatrix PearsonCorr3.tab
echo "1z done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_PearsonCorr4.png --outFileCorMatrix PearsonCorr4.tab
echo "2z done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_SpearmanCorr3.png --outFileCorMatrix SpearmanCorr3.tab
echo "3z done"
/home/wangjiachen/anaconda3/envs/deeptools2/bin/plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_SpearmanCorr4.png --outFileCorMatrix SpearmanCorr4.tab
echo "4z done"
echo "plot correlation done"
mkdir Correlationplot ; mv *.tab Correlationplot ; mv *.png Correlationplot ; mv *.npz Correlationplot
echo "All done"

