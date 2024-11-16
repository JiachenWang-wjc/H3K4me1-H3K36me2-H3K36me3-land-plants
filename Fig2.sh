#Figure 2A (Take rice for example)
computeMatrix scale-regions -S \
WT_H3K4me1.bigwig \
WT_RNA.bigwig \
-R Genes.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -p 30 -out H3K4me1.mat.gz
plotHeatmap -m H3K4me1.mat.gz -out H3K4me1_heatmap.pdf --sortRegions descend --sortUsingSamples 1 --colorList "white,#FF0000" --zMax 50


#Figure 2B
plotProfile -m H3K4me1.mat.gz   -out H3K4me1_pattern.pdf  --perGroup --colors black red --samplesLabel NIP --legendLocation upper-right --yMin 0 --yMax 50 --plotHeight 11 --plotWidth 11


#Figure 2C
FPKM <- read.delim("FPKM_NIP_all.xls")
FPKM <- subset(FPKM,FPKM$NIP_mean >= 1)
H3K4me1_enriched_genes <- read.delim("H3K4me1peakRelatedgene_annoOvgene_genelist.xls")
row.names(H3K4me1_enriched_genes) <- H3K4me1_enriched_genes[,1]
target <- intersect(row.names(H3K4me1_enriched_genes),FPKM[,1])
FPKM <- subset(FPKM, is.element(FPKM[,1], target) == T)
FPKM <- FPKM[order(FPKM[,6],decreasing = F),]
library(dplyr)
FPKM <- FPKM %>% mutate(group = ntile(FPKM$NIP_mean, 100))
input <- read.delim("input.mat", header = F)
raw <- read.delim("H3K36me3.mat", header = F)
row.names(input) <- input[,4]
row.names(raw) <- raw[,4]
H3K4me1 <- raw[,7:606] - input[,7:606]
H3K4me1 <- H3K4me1[, 151:450]
H3K4me1$mean <- rowMeans(H3K4me1)
H3K4me1 <- subset(H3K4me1, is.element(row.names(H3K4me1), target) == T)
FPKM <- FPKM[order(FPKM[,1], decreasing = F),]
H3K4me1 <- H3K4me1[order(row.names(H3K4me1), decreasing = F),]
data <- data.frame(FPKM[,16], H3K4me1[,301])
row.names(data) <- row.names(H3K4me1)
names(data) <- c("FPKM_group", "methylation")
data1 <- data %>% group_by(FPKM_group) %>% summarise_all(mean)
data1 <- data1[order(data1$methylation, decreasing = F),]
library(dplyr)
data1 <- data1 %>% mutate(methylation_group = ntile(data1$methylation, 100))
data1 <- data1[, c(1,3)]
data1 <- data1[order(data1$FPKM_group, decreasing = F),]
cor.test(data1$FPKM_group, data1$methylation_group, alternative = "two.sided", method = "pearson")
cor.test(data1$FPKM_group, data1$methylation_group, alternative = "two.sided", method = "spearman")
library("ggplot2")
library("showtext")
showtext_auto(enable = TRUE)
font_add('TNR', 'times.ttf')
ggplot(data1, aes(x = FPKM_group, y = methylation_group)) +
  geom_point() + theme(legend.position='none') +
  theme(panel.grid=element_blank(), panel.background=element_rect(color="black",size = 1, fill = "white")) +
  labs(x = "Rank of gene expression levels",y = "Rank of methylation levels") +
  theme(axis.title.x= element_text(size=20, family="TNR",color="black", vjust=0.5, hjust=0.5)) +
  theme(axis.text.x =  element_text(size=18, family="TNR",color="black",vjust=0.5, hjust=0.5)) +
  theme(axis.title.y= element_text(size=20, family="TNR",color="black", vjust=0.5, hjust=0.5)) +
  theme(axis.text.y =  element_text(size=18, family="TNR",color="black",vjust=0.5, hjust=0.5)) 











