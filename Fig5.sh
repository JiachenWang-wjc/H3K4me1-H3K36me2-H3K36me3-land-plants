#Figure 5G
#Find DEGs
DEGfind=function(controlname,treatname,fFoldcp,fFoldcq,pPalue,pPadj) 
{
dirnow <- getwd()
rawcount <- data.frame(row.names=sp1[,1], control_1=sp1[,2],control_2=sp2[,2],control_3=sp3[,2],treat_1=sp4[,2],treat_2=sp5[,2],treat_3=sp6[,2])
condition <- factor(rep(c("control","treat"), each = 3))
suppressMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(rawcount, DataFrame(condition), design = ~condition)
dds_filter <- dds[rowSums(counts(dds)) >= 1,]
dds_DES <- DESeq(dds_filter)
res <- results(dds_DES, contrast = c("condition", "treat", "control"))
fFoldcp <- as.numeric(fFoldcp);pPalue <- as.numeric(pPalue)
dir.create(paste(treatname,"vs",controlname,"FC",fFoldcp,"p",pPalue,sep="_"));dirname <- paste(treatname,"vs",controlname,"FC",fFoldcp,"p",pPalue,sep="_")
setwd(dirname)
res$FC <- 2^res$log2FoldChange
resa <- res[res$FC<1,];resa$FCx <- -(1/resa$FC);resb <- res[res$FC>=1,];resb$FCx <- resb$FC;resc <- rbind(resa,resb)
res1 <- resc[order(resc$pvalue),];res2<-res1[,-3];res3<-res2[,-3]
res1b <- resc[order(resc$padj),];res2b<-res1b[,-3];res3b<-res2b[,-3]
diff <- subset(res3, abs(res3$log2FoldChange) >= log2(fFoldcp) & res3$pvalue <= pPalue)
up <- subset(res3, res3$log2FoldChange >= log2(fFoldcp) & res3$pvalue <= pPalue)
down <- subset(res3, res3$log2FoldChange <= -log2(fFoldcp) & res3$pvalue <= pPalue)
Allgenes<-row.names(res3);diffgenes<-row.names(diff);upgenes<-row.names(up);downgenes<-row.names(down)
res0 <- cbind(Allgenes,res3);diff <- cbind(diffgenes,diff);up <- cbind(upgenes,up);down <- cbind(downgenes,down)
lowexpress <- setdiff(sp1[,1],res0[,1])
lowexpress <- as.data.frame(lowexpress)
lowexpress$Allgenes <- lowexpress[,1]
lowexpress$baseMean <- 0;lowexpress$log2FoldChange <- 0;lowexpress$pvalue <- NA;lowexpress$padj <- NA;lowexpress$FC <- 1;lowexpress$FCx <- -1
lowexpress <- lowexpress[,-1]
resppp <- as.data.frame(res0[,1]);resppp$Allgenes <- resppp[,1]
resppp$baseMean <- res0$baseMean;resppp$log2FoldChange <- res0$log2FoldChange;resppp$pvalue <- res0$pvalue;resppp$padj <- res0$padj;resppp$FC <- res0$FC;resppp$FCx <- res0$FCx
resppp <- resppp[,-1]
res0A <- rbind(resppp,lowexpress)
write.table(res0A,paste(treatname,fFoldcp,pPalue,"vs",controlname,"p.xls",sep="_"),sep="\t",row.names=F)
write.table(diff,paste(treatname,fFoldcp,pPalue,"diff",controlname,"p.xls",sep="_"),sep="\t",row.names=F)
write.table(up,paste(treatname,fFoldcp,pPalue,"up",controlname,"p.xls",sep="_"),sep="\t",row.names=F)
write.table(down,paste(treatname,fFoldcp,pPalue,"down",controlname,"p.xls",sep="_"),sep="\t",row.names=F)
setwd(dirnow)
pPadj <- as.numeric(pPadj);fFoldcq <- as.numeric(fFoldcq)
dir.create(paste(treatname,"vs",controlname,"FC",fFoldcq,"padj",pPadj,sep="_"))
dirname1 <- paste(treatname,"vs",controlname,"FC",fFoldcq,"padj",pPadj,sep="_");setwd(dirname1)
diff1 <- subset(res3b, abs(res3b$log2FoldChange) >= log2(fFoldcq) & res3b$padj <= pPadj)
up1 <- subset(res3b, res3b$log2FoldChange >= log2(fFoldcq) & res3b$padj <= pPadj)
down1 <- subset(res3b, res3b$log2FoldChange <= -log2(fFoldcq) & res3b$padj <= pPadj)
Allgenes1<-row.names(res3b);diffgenes1<-row.names(diff1);upgenes1<-row.names(up1);downgenes1<-row.names(down1)
res00 <- cbind(Allgenes1,res3b);diff1 <- cbind(diffgenes1,diff1);up1 <- cbind(upgenes1,up1);down1 <- cbind(downgenes1,down1)
lowexpressA <- setdiff(sp1[,1],res00[,1])
lowexpressA <- as.data.frame(lowexpressA)
lowexpressA$Allgenes1 <- lowexpress[,1]
lowexpressA$baseMean <- 0;lowexpressA$log2FoldChange <- 0;lowexpressA$pvalue <- NA;lowexpressA$padj <- NA;lowexpressA$FC <- 1;lowexpressA$FCx <- -1
lowexpressA <- lowexpressA[,-1]
resbbb <- as.data.frame(res00[,1]);resbbb$Allgenes1 <- resbbb[,1]
resbbb$baseMean <- res00$baseMean;resbbb$log2FoldChange <- res00$log2FoldChange;resbbb$pvalue <- res00$pvalue;resbbb$padj <- res00$padj;resbbb$FC <- res00$FC;resbbb$FCx <- res00$FCx
resbbb <- resbbb[,-1]
res0B <- rbind(resbbb,lowexpressA)
write.table(res0B,paste(treatname,fFoldcq,pPadj,"vs",controlname,"padj.xls",sep="_"),sep="\t",row.names=F)
write.table(diff1,paste(treatname,fFoldcq,pPadj,"diff",controlname,"padj.xls",sep="_"),sep="\t",row.names=F)
write.table(up1,paste(treatname,fFoldcq,pPadj,"up",controlname,"padj.xls",sep="_"),sep="\t",row.names=F)
write.table(down1,paste(treatname,fFoldcq,pPadj,"down",controlname,"padj.xls",sep="_"),sep="\t",row.names=F)
setwd(dirnow);nrow(diff)
}
sp1  <- read.delim("merge_NiP_14days_rep1_rmCMUnSy_q20_s.rawcount")
sp2  <- read.delim("merge_NiP_14days_rep2_rmCMUnSy_q20_s.rawcount")
sp3  <- read.delim("merge_NiP_14days_rep3_rmCMUnSy_q20_s.rawcount")
sp4  <- read.delim("merge_ehd3_14days_rep1_rmCMUnSy_q20_s.rawcount")
sp5  <- read.delim("merge_ehd3_14days_rep2_rmCMUnSy_q20_s.rawcount")
sp6  <- read.delim("merge_ehd3_14days_rep3_rmCMUnSy_q20_s.rawcount")
DEGfind("NiP","ehd3","1.5","1.5","0.05","0.05")
sp4  <- read.delim("merge_sdg724_14days_rep1_rmCMUnSy_q20_s.rawcount")
sp5  <- read.delim("merge_sdg724_14days_rep2_rmCMUnSy_q20_s.rawcount")
sp6  <- read.delim("merge_sdg724_14days_rep3_rmCMUnSy_q20_s.rawcount")
DEGfind("NiP","sdg724","1.5","1.5","0.05","0.05")
sp4  <- read.delim("merge_ehd3sdg724_14days_rep1_rmCMUnSy_q20_s.rawcount")
sp5  <- read.delim("merge_ehd3sdg724_14days_rep2_rmCMUnSy_q20_s.rawcount")
sp6  <- read.delim("merge_ehd3sdg724_14days_rep3_rmCMUnSy_q20_s.rawcount")
DEGfind("NiP","ehd3sdg724","1.5","1.5","0.05","0.05")
#Heatmap (Adjust the column information and rename the output from the previous step.)
Sample1_all <- read.delim("ehd3 vs WT padj0.05 log2(1.5).xls")
Sample2_all <- read.delim("sdg724 vs WT padj0.05 log2(1.5).xls")
Sample3_all <- read.delim("ehd3sdg724 vs WT padj0.05 log2(1.5).xls")
Sample1_diff <- read.delim("ehd3 diff WT padj0.05 log2(1.5).xls")
Sample2_diff <- read.delim("sdg724 diff WT padj0.05 log2(1.5).xls")
Sample3_diff <- read.delim("ehd3sdg724 diff WT padj0.05 log2(1.5).xls")
genelist <- unique(c(row.names(Sample1_diff),row.names(Sample2_diff),row.names(Sample3_diff)))
length(genelist)
Sample1_tar <-  subset(Sample1_all,is.element(row.names(Sample1_all),genelist)==T)  
Sample2_tar <-  subset(Sample2_all,is.element(row.names(Sample2_all),genelist)==T)  
Sample3_tar <-  subset(Sample3_all,is.element(row.names(Sample3_all),genelist)==T)  
countdata <- data.frame(ehd3=Sample1_tar$log2FoldChange,sdg724=Sample2_tar$log2FoldChange,ehd3sdg724=Sample3_tar$log2FoldChange,
                        row.names=row.names(Sample1_tar))
write.table(countdata,"heatmap_countdata.xls",sep="\t")
library(ComplexHeatmap)
library(circlize)
countdata <- read.delim("heatmap_countdata.xls")
heatdata <- countdata
k=kmeans(heatdata,2)
heatdata$k=k$cluster
htframe_sort=heatdata[order(heatdata$k),]
Heatmap(htframe_sort[,1:(ncol(htframe_sort)-1)],col=colorRamp2(c(-log2(4), 0, log2(4)), c("#4779BD", "white", "red")),show_row_names=FALSE,cluster_rows=FALSE,cluster_columns=FALSE, width = unit(5, "cm"))
write.table(htframe_sort,"htframe_sort_1.xls",sep="\t")
htframe_sort2 <- rbind(subset(htframe_sort,htframe_sort[,6]==2),subset(htframe_sort,htframe_sort[,6]==1))
write.table(htframe_sort2,"htframe_sort_2.xls",sep="\t")
h1<-Heatmap(htframe_sort2[,1:(ncol(htframe_sort)-1)],col=colorRamp2(c(-log2(4), 0, log2(4)), c("#4779BD", "white", "red")),show_row_names=FALSE,cluster_rows=FALSE,cluster_columns=FALSE, width = unit(3, "cm"))
ha1_row = rowAnnotation(df = data.frame(type2 =htframe_sort2$k),
         col = list(type2 = c("1" = "blue","2"="red")), width = unit(0.3, "cm"))
draw(ha1_row+h1)


#Figure 5H
volcano_correlation <- function(Path,file1,file2,log2FC,p_value,positive_value,negative_value,pdfname,xlab,ylab,xrange,yrange) {
  setwd(Path)
  library(ggplot2)
  gene1 <- read.delim(file1)
  gene2 <- read.delim(file2)
  gene1 <- gene1 [order(row.names(gene1),decreasing=F),]
  gene2 <- gene2 [order(row.names(gene2),decreasing=F),]
  all <- data.frame(gene1 [], gene2 [])
  all_gene1_up   <- subset(all, all$log2FoldChange >= log2FC & all$padj <= p_value )
  all_gene1_down <- subset(all, all$log2FoldChange <= -log2FC & all$padj <= p_value )
  all_gene2_up   <- subset(all, all$log2FoldChange.1 >= log2FC & all$padj.1 <= p_value )
  all_gene2_down <- subset(all, all$log2FoldChange.1 <= -log2FC & all$padj.1 <= p_value )
  print(c(nrow(all_gene1_up),nrow(all_gene1_down),nrow(all_gene2_up),nrow(all_gene2_down)))
  allbind <- rbind(all_gene1_up,all_gene1_down,all_gene2_up,all_gene2_down)
  all <- unique(allbind)
  gene_up <- subset(all, all$log2FoldChange >= log2FC & all$padj <= p_value & all$log2FoldChange.1 >= log2FC & all$padj.1 <= p_value )
  gene_down <- subset(all, all$log2FoldChange <= -log2FC & all$padj <= p_value & all$log2FoldChange.1 <= -log2FC & all$padj.1 <= p_value )
  print(c(nrow(gene_up), nrow(gene_down)))
  gene_up_location <- Reduce(intersect,list(v1=which(all$log2FoldChange >= log2FC), v2=which(all$padj <= p_value), v3=which(all$log2FoldChange.1 >= log2FC), v4= which(all$padj.1 <= p_value) )) 
  gene_down_location <- Reduce(intersect,list(v1=which(all$log2FoldChange <= -log2FC), v2=which(all$padj <= p_value), v3=which(all$log2FoldChange.1 <= -log2FC), v4= which(all$padj.1 <= p_value) )) 
  sig <- rep("no", times = nrow(all))
  sig[gene_up_location] <- "up"
  sig[gene_down_location] <- "down"  
  sig <- factor(sig,levels=c("up","down","no")) 
  all$type <- sig
  xline <- c(-log2FC,log2FC)
  yline <-  c(-log2FC,log2FC)
  cor <- cor.test(all$log2FoldChange, all$log2FoldChange.1, alternative = "two.sided", method = "pearson")
  print(cor)
  all$log2FoldChange[all$log2FoldChange > positive_value] <- positive_value
  all$log2FoldChange[all$log2FoldChange < negative_value] <- negative_value
  all$log2FoldChange.1[all$log2FoldChange.1 > positive_value] <- positive_value
  all$log2FoldChange.1[all$log2FoldChange.1 < negative_value] <- negative_value
  p <- ggplot(all, aes(x = log2FoldChange,y = log2FoldChange.1, colour = type)) +
    geom_point() +
    theme(legend.position='none') +
    theme(panel.grid=element_blank(), panel.background=element_rect(color="black",size = 1, fill = "white")) +
    scale_discrete_manual(values = c("red","green","black"), aesthetics = 'colour') +
    geom_hline(yintercept = yline, lty = 2, size = I(0.5), colour = "grey") + 
    geom_vline(xintercept = xline, lty = 2,size = I(0.5), colour = "grey") +
    labs(x = xlab,y = ylab) +
    theme(axis.title.x= element_text(size=20, color="black", vjust=0.5, hjust=0.5)) +
    theme(axis.text.x =  element_text(size=18, color="black",vjust=0.5, hjust=0.5)) +
    theme(axis.title.y= element_text(size=20, color="black", vjust=0.5, hjust=0.5)) +
    theme(axis.text.y =  element_text(size=18, color="black",vjust=0.5, hjust=0.5)) + 
    scale_x_continuous(limits = c( -xrange,yrange)) +
    scale_y_continuous(limits = c( -yrange,yrange)) 
  ggsave(p, file=pdfname,width = 5,height = 5) 
}
volcano_correlation("C:/Users/Desktop","ehd3 vs WT padj0.05 log2(1.5).xls", "sdg724 vs WT padj0.05 log2(1.5).xls",log2(1.5),0.05,10,-10,"volcano_correlation ehd3 vs sdg724.pdf","log2FC  ehd3/WT","log2FC  sdg724/WT",10,10)
#cor=0.7437783 p-value < 2.2e-16
volcano_correlation("C:/Users/Desktop","ehd3 vs WT padj0.05 log2(1.5).xls", "ehd3sdg724 vs WT padj0.05 log2(1.5).xls",log2(1.5),0.05,10,-10,"volcano_correlation ehd3 vs ehd3sdg724.pdf","log2FC  ehd3/WT","log2FC  ehd3sdg724/WT",10,10)
#cor=0.8506074 p-value < 2.2e-16
volcano_correlation("C:/Users/Desktop","sdg724 vs WT padj0.05 log2(1.5).xls", "ehd3sdg724 vs WT padj0.05 log2(1.5).xls",log2(1.5),0.05,10,-10,"volcano_correlation sdg724 vs ehd3sdg724.pdf","log2FC  sdg724/WT","log2FC  ehd3sdg724/WT",10,10)
#cor=0.8288915 p-value < 2.2e-16


#Figure 5I
##Linux shell
#call peak (Take one of these samples for example.) with Input control using SICER or MACS3.

#SICER
sh SICER.sh /Path NiP-K36me1.bed NiP-Input.bed /OutputPath msu7  1 200 200 0.9 600 0.05
#Results with FDR<0.05 were used as thresholds. 
#MACS
macs3 callpeak -t IP.bam -c Input.bam -n IP_q0.05_broad --outdir IP_q0.05_broad -f BAMPE -g 3.7e8 --broad --broad-cutoff 0.05 -q 0.05  --nomodel  -B
#Name the peaks result as a new BED file.

#The intersection of different biological repeats was used to obtain stable peaks.
##R
Peakover <- function(Afilename,Bfilename,peak1name,peak2name,mainname){
 #packages
 suppressMessages(library("ChIPpeakAnno"))
 suppressMessages(library("GenomicRanges"))
 #Afilename-filename of Afile; Bfilename-filename of Bfile; peak1name-name of peaks in Afile; peak2name-name of peaks in Bfile; mainname-prefiex of name in outputfile
 Afile <- read.delim(Afilename,header=F)
 Bfile <- read.delim(Bfilename,header=F)
 Agranges <- GRanges(seqnames=Afile[,1],ranges=IRanges(start=Afile[,2],end=Afile[,3],names=paste("A",rep(1:nrow(Afile)),sep="")))
 Bgranges <- GRanges(seqnames=Bfile[,1],ranges=IRanges(start=Bfile[,2],end=Bfile[,3],names=paste("B",rep(1:nrow(Bfile)),sep="")))
 result <<- findOverlapsOfPeaks(Agranges,Bgranges,connectedPeaks="merge")
 pdf(paste(mainname,".pdf",sep=""))
 makeVennDiagram(result,NameOfPeaks=c(peak1name,peak2name),height=3000,width=3000,col="transparent",fill=c("red","green"),alpha=c(0.5,0.5),main=mainname)
 dev.off()
 write.table(data.frame(result$mergedPeaks)[,1:3],paste(mainname,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F) 
 }

Peakover("NiP_H3K4me1_rep1_q0.05_broad_peaks.broadPeak","NiP_H3K4me1_rep2_q0.05_broad_peaks.broadPeak","NiP_H3K4me1_rep1","NiP_H3K4me1_rep2","NiP_H3K4me1_1v2")
Peakover("NiP_H3K4me1_rep3_q0.05_broad_peaks.broadPeak","NiP_H3K4me1_1v2.bed","NiP_H3K4me1_rep3","temp","NiP_H3K4me1_merge")

#Diff peak
##Linux shell
manorm  --p1 mutant_K36me1_peak.bed --p2 NiP_K36me1_peak.bed --pf bed --r1 mutant-K36me1.bam --r2 NiP-K36me1.bam \
--rf bam --pe -w 2000 -m 0.26 -p 0.05 --name1 mutant --name2 NiP -o K36me1_mutant_diffpeak
#Heatmap from MAnorm result
##R
library(ggplot2);library(dplyr);library(viridis);library(ggpointdensity);library(cowplot)
#H3K4me1-ehd3
H3K4me1_ehd3_rep1 <- read.delim("NiP_H3K4me1_rep1_vs_ehd3_H3K4me1_rep1_all_MAvalues.xls", header = T)
H3K4me1_ehd3_rep2 <- read.delim("NiP_H3K4me1_rep2_vs_ehd3_H3K4me1_rep2_all_MAvalues.xls", header = T)
H3K4me1_ehd3_rep3 <- read.delim("NiP_H3K4me1_rep3_vs_ehd3_H3K4me1_rep3_all_MAvalues.xls", header = T)
H3K4me1_ehd3 <- (H3K4me1_ehd3_rep1[, 9:10] + H3K4me1_ehd3_rep2[, 9:10] + H3K4me1_ehd3_rep3[, 9:10])/300
p1<- ggplot(data = H3K4me1_ehd3, mapping = aes(x = H3K4me1_ehd3[,1], y = H3K4me1_ehd3[,2]))+theme_classic()+geom_pointdensity(size = 0.5)+scale_color_viridis(option = "G")
p1 <- p1+xlim(0,3.5)+ylim(0,3.5)+xlab("NiP")+ylab("ehd3")+ ggtitle("H3K4me1_Normalized_Level")+geom_abline(intercept = 0, slope = 1,linetype = "dashed", size = 1)
p1