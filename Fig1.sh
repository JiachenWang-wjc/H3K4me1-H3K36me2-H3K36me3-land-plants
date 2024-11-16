#Figure 1A (Take the heatmap of H3K36me2 for example)
computeMatrix scale-regions -S WT_H3K36me2_merge_for_sdg724.bw WT_Input_merge_for_sdg724.bw sdg724_H3K36me2_merge.bw sdg724_Input_merge.bw \
WT_H3K36me2_merge_for_sdg725.bw WT_Input_merge_for_sdg725.bw sdg725_H3K36me2_merge.bw sdg725_Input_merge.bw \
WT_H3K36me2_merge_for_sdg701.bw WT_Input_merge_for_sdg701.bw sdg701_H3K36me2_merge.bw sdg701_Input_merge.bw \
WT_H3K36me2_merge_for_osino80.bw WT_Input_merge_for_osino80.bw osino80_H3K36me2_merge.bw osino80_Input_merge.bw \
-R NiP_H3K36me2_enriched_genes.bed -b 3000 -a 3000 -m 6000 --binSize 10 --sortRegions keep --missingDataAsZero -p 45 -out H3K36me2.mat.gz
R
allm <- read.delim("H3K36me2.mat", header = F, skip=1);row.names(allm) <- allm[,4];nrow(allm);ncol(allm)
WT1_H3K36me2 <- allm[, 7:1206];WT1_ipt <- allm[, 1207:2406];sdg724_H3K36me2 <- allm[, 2407:3606];sdg724_ipt <- allm[, 3607:4806]
WT2_H3K36me2 <- allm[, 4807:6006];WT2_ipt <- allm[, 6007:7206];sdg725_H3K36me2 <- allm[, 7207:8406];sdg725_ipt <- allm[, 8407:9606]
WT3_H3K36me2 <- allm[, 9607:10806];WT3_ipt <- allm[, 10807:12006];sdg701_H3K36me2 <- allm[, 12007:13206];sdg701_ipt <- allm[, 13207:14406]
WT4_H3K36me2 <- allm[, 14407:15606];WT4_ipt <- allm[, 15607:16806];osino80_H3K36me2 <- allm[, 16807:18006];osino80_ipt <- allm[, 18007:19206]
xy1<-sum(WT1_H3K36me2);print(xy1)
xy2<-sum(WT1_ipt);print(xy2)
xy3<-sum(sdg724_H3K36me2);print(xy3)
xy4<-sum(sdg724_ipt);print(xy4)
xy5<-sum(WT2_H3K36me2);print(xy5)
xy6<-sum(WT2_ipt);print(xy6)
xy7<-sum(sdg725_H3K36me2);print(xy7)
xy8<-sum(sdg725_ipt);print(xy8)
xy9<-sum(WT3_H3K36me2);print(xy9)
xy10<-sum(WT3_ipt);print(xy10)
xy11<-sum(sdg701_H3K36me2);print(xy11)
xy12<-sum(sdg701_ipt);print(xy12)
xy13<-sum(WT4_H3K36me2);print(xy13)
xy14<-sum(WT4_ipt);print(xy14)
xy15<-sum(osino80_H3K36me2);print(xy15)
xy16<-sum(osino80_ipt);print(xy16)
WT1_H3K36me2 <- 10^9 * (WT1_H3K36me2)/xy1
WT1_ipt <- 10^9 * (WT1_ipt)/xy2
sdg724_H3K36me2 <- 10^9 * (sdg724_H3K36me2)/xy3
sdg724_ipt <- 10^9 * (sdg724_ipt)/xy4
WT2_H3K36me2 <- 10^9 * (WT2_H3K36me2)/xy5
WT2_ipt <- 10^9 * (WT2_ipt)/xy6
sdg725_H3K36me2 <- 10^9 * (sdg725_H3K36me2)/xy7
sdg725_ipt <- 10^9 * (sdg725_ipt)/xy8
WT3_H3K36me2 <- 10^9 * (WT3_H3K36me2)/xy9
WT3_ipt <- 10^9 * (WT3_ipt)/xy10
sdg701_H3K36me2<- 10^9 * (sdg701_H3K36me2)/xy11
sdg701_ipt <- 10^9 * (sdg701_ipt)/xy12
WT4_H3K36me2 <- 10^9 * (WT4_H3K36me2)/xy13
WT4_ipt <- 10^9 * (WT4_ipt)/xy14
osino80_H3K36me2 <- 10^9 * (osino80_H3K36me2)/xy15
osino80_ipt <- 10^9 * (osino80_ipt)/xy16
geneBed <- allm[,1:6]
WT1_H3K36me2 <- WT1_H3K36me2 - WT1_ipt
sdg724_H3K36me2 <- sdg724_H3K36me2 - sdg724_ipt
WT2_H3K36me2 <- WT2_H3K36me2 - WT2_ipt
sdg725_H3K36me2 <- sdg725_H3K36me2 - sdg725_ipt
WT3_H3K36me2 <- WT3_H3K36me2 - WT3_ipt
sdg701_H3K36me2 <- sdg701_H3K36me2 - sdg701_ipt
WT4_H3K36me2 <- WT4_H3K36me2 - WT4_ipt
osino80_H3K36me2 <- osino80_H3K36me2 - osino80_ipt
WT1_H3K36me2[WT1_H3K36me2 < 0] <- 0
sdg724_H3K36me2[sdg724_H3K36me2 < 0] <- 0
WT2_H3K36me2[WT2_H3K36me2 < 0] <- 0
sdg725_H3K36me2[sdg725_H3K36me2 < 0] <- 0
WT3_H3K36me2[WT3_H3K36me2 < 0] <- 0
sdg701_H3K36me2[sdg701_H3K36me2 < 0] <- 0
WT4_H3K36me2[WT4_H3K36me2 < 0] <- 0
osino80_H3K36me2[osino80_H3K36me2 < 0] <- 0
change_sdg724 <- log2((sdg724_H3K36me2+1) / (WT1_H3K36me2+1))
change_sdg725 <- log2((sdg725_H3K36me2+1) / (WT2_H3K36me2+1))
change_sdg701 <- log2((sdg701_H3K36me2+1) / (WT3_H3K36me2+1))
change_osino80<- log2((osino80_H3K36me2+1) / (WT4_H3K36me2+1))
df <- cbind(geneBed,WT1_H3K36me2,change_sdg724,change_sdg725,change_sdg701,change_osino80)
write.table(df,"BPM_Change_H3K36me2.mat",sep="\t",row.names=F,col.names=F,quote=F);nrow(df);ncol(df)
#Linux-shell
sed -i '1 i @{"upstream":[3000,3000,3000,3000,3000],"downstream":[3000,3000,3000,3000,3000],"body":[6000,6000,6000,6000,6000],"bin size":[10,10,10,10,10],"ref point":[null,null,null,null,null],"verbose":false,"bin avg type":"mean","missing data as zero":true,"min threshold":null,"max threshold":null,"scale":1,"skip zeros":false,"nan after end":false,"proc number":4,"sort regions":"keep","sort using":"mean","unscaled 5 prime":[0,0,0,0,0],"unscaled 3 prime":[0,0,0,0,0],"group_labels":["n21085"],"group_boundaries":[0,21085],"sample_labels":["NiP_H3K36me2","change_sdg724","change_sdg725","change_sdg701","change_osino80"],"sample_boundaries":[0,1200,2400,3600,4800,6000]}' BPM_Change_H3K36me2.mat
gzip -1 BPM_Change_H3K36me2.mat
plotHeatmap -m BPM_Change_H3K36me2.mat.gz -out BPM_Change_H3K36me2.pdf --outFileSortedRegions BPM_Change_H3K36me2.bed --sortUsingSamples 1 --sortRegions descend  --colorList "white,#409E76" "#104E8B,white,#CD3700" "#104E8B,white,#CD3700" "#104E8B,white,#CD3700" "#104E8B,white,#CD3700" --samplesLabel  NiP_H3K36me2 change_sdg724 change_sdg725 change_sdg701 change_osino80 --zMin 0 -4 -4 -4 -4  --zMax 50 4 4 4 4  --yMin 0 --yMax 50


#Figure 1B
#
#This example takes ehd3 as a case study. If you want to add multiple samples, modify the corresponding code according to the number of samples and their names.
#
#1.Number the exons of each gene.
##R
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("all_nounsy.gtf", format = "gtf")
#Find the exons corresponding to their genes.
geneexonlists<-exonsBy(txdb, "gene", use.names=F)
write.table(geneexonlists,"Genes_and_exons.xls",sep="\t",row.names=F)
#Use the following Microsoft EXCEL function to sort the exons of each gene and add a new name column.
# "+" strand:COUNTIF(B$2:B2,B2)
# "-" strand:COUNTIF(B2:B$127410, B2);####127410 is the total number of exons on the "-" strand (plus one).
#They were stored separately as bed files (a for "+" ; b for "-").

##Linux shell
#Sort the bed files by chromosome and then by start position.
sort -k1,1 -k2,2n rice_exon_nomerge_a.bed > rice_exon_nomerge_a_s.bed
sort -k1,1 -k2,2n rice_exon_nomerge_b.bed > rice_exon_nomerge_b_s.bed

#2.Find exons that overlap on the genome and merge them (Find all exons belonging to the same gene for each gene).
bedtools merge -s -S + -c 1,2,3,4,5,6 -o  distinct,collapse,collapse,collapse,collapse,count -i rice_exon_nomerge_a_s.bed > rice_exon_nomerge_a_ss.bed
bedtools merge -s -S - -c 1,2,3,4,5,6 -o  distinct,collapse,collapse,collapse,collapse,count -i rice_exon_nomerge_b_s.bed > rice_exon_nomerge_b_ss.bed
#Using the EXCEL function above, rename the non-overlapping exons.
#Combine the two files and get a full list for verification and visualization.
cat rice_exon_merge_a.bed rice_exon_merge_b.bed > rice_exon_merge.bed
sort -k1,1 -k2,2n rice_exon_merge.bed > rice_exon_merge_s.bed

#3.Count the reads number on each exon.
bedtools multicov -D -bams WT_merge.bam ehd3_merge.bam -bed rice_exon_merge_a.bed > count_a.txt
bedtools multicov -D -bams WT_merge.bam ehd3_merge.bam -bed rice_exon_merge_b.bed > count_b.txt
#Biological replicates could be merged in multiple steps, for example, by using bedtools to count reads for each library individually. After testing, the choice of different merging orders makes very little difference.
cat count_a.txt count_b.txt > count.txt


#4.Calculation of exon-based TPM (transcripts per million).
#Remove the location information and prepare to calculate exon-based TPM (Store as a new bed file).
##R
#Calculation of exon-based TPM.
ecount  <- read.delim("count_new.txt", header=F);nrow(ecount);ncol(ecount)
row.names(ecount) <- ecount[, 1]
rate_WT1 <-ecount[, 5]/ecount[, 3]
rate_ehd3 <-ecount[, 6]/ecount[, 3]
WT1_total_rate <- sum(rate_WT1)
ehd3_total_rate <- sum(rate_ehd3)
TPM_WT1 <- rate_WT1*(10^6)/WT1_total_rate
TPM_ehd3 <- rate_ehd3*(10^6)/ehd3_total_rate
TPM <- data.frame(TPM_WT1,TPM_ehd3)
Allexons<-row.names(ecount)
Allgenes<-ecount[, 2]
eTPM <- cbind(Allexons,Allgenes,TPM)
write.table(eTPM,"eTPM_all.txt",sep="\t",row.names=F, quote=FALSE)

#5.Perform eFC (exon-based Fold of chang of TPM) calculations.
#H3K4me1 enriched genes (n=17228) were selected from the list above.
##Linux shell
awk -F'\t' 'NR==FNR {keys[$1]; next} $2 in keys' H3K4me1_enriched_genes_n17228.txt eTPM_all.txt > eTPM_all_n17228.txt
#R
#eFC=log2((eTPM_mutant+ε)/(eTPM_WT+ε));ε=1
eTPM_n17228  <- read.delim("eTPM_all_n17228.txt", header=F)
nrow(eTPM_n17228);ncol(eTPM_n17228)
row.names(eTPM_n17228) <- eTPM_n17228[, 1]
eFC_ehd3 <-log2((eTPM_n17228[, 4]+1)/(eTPM_n17228[, 3]+1))
eFC <- data.frame(eFC_ehd3)
n17228exon_names<-row.names(eTPM_n17228)
n17228gene_names<-eTPM_n17228[, 2]
eFC_new <- cbind(n17228exon_names,n17228gene_names,eFC)
write.table(eFC_new,"eFC_n17228.txt",sep="\t",row.names=F, quote=FALSE);nrow(eFC_new);ncol(eFC_new)
##Linux shell
#Extract the first, second, third, fourth, fifth, tenth, fifteenth, and twentieth exons and their eFCs.
awk -F'\t' '$1 ~ /\.1$/' eFC_n17228.txt >eFC_N1_1st.xls
awk -F'\t' '$1 ~ /\.2$/' eFC_n17228.txt >eFC_N1_2nd.xls
awk -F'\t' '$1 ~ /\.3$/' eFC_n17228.txt >eFC_N1_3rd.xls
awk -F'\t' '$1 ~ /\.4$/' eFC_n17228.txt >eFC_N1_4th.xls
awk -F'\t' '$1 ~ /\.5$/' eFC_n17228.txt >eFC_N1_5th.xls
awk -F'\t' '$1 ~ /\.10$/' eFC_n17228.txt >eFC_N1_10th.xls
awk -F'\t' '$1 ~ /\.15$/' eFC_n17228.txt >eFC_N1_15th.xls
awk -F'\t' '$1 ~ /\.20$/' eFC_n17228.txt >eFC_N1_20th.xls
#Paste the output information into a new tab-separated txt file and add the exon serial number in the first line.

#6.Box-plot
library(tidyverse)
datag <- read.delim('eFC_boxplot.txt',header = T)
datag <- gather(datag) %>% dplyr::rename('Sample'='key','log2eFC'='value')
clean_datag <- na.omit(datag)
ggplot(clean_datag,aes(`Sample`,`log2eFC`))+stat_boxplot(geom="errorbar",width=0.5,size=0.5)+geom_boxplot(width=0.75,cex=0.5,outlier.shape=NA,aes(fill=Sample))+ylim(-0.5,0.5)+
geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8)+
scale_fill_manual(values = c("#76EEC6","#76EEC6","#76EEC6","#76EEC6","#76EEC6","#76EEC6","#76EEC6","#76EEC6"))+
theme_bw()+theme(panel.background=element_rect(colour="black",size=1,fill="white"),axis.line = element_line(colour = "black",size=0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(colour="black",face="bold",size=14,angle = 45,vjust = 0.5,hjust = 0.5), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=14), 
        axis.title.x=element_text(size=15,face="bold"),
        axis.title.y=element_text(size=15,face="bold")) +
  labs(x = "Group", y = "log2eFC")     
#Save the PDF file as 50*10inches.