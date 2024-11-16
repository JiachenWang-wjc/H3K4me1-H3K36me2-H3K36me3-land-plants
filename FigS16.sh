#差异基因火山图2
#volcanoplot ehd3 vs WT padj0.05 log2(1.5)
setwd(dir = "C:/Users/HP/Desktop")
gene <- read.delim("E:/复旦/实验室/课题/8_炳师兄水稻/new/RNAseq/14days/差异基因/合并计算/padj0.05 log2(1.5)/ehd3 vs WT padj0.05 log2(1.5).xls")
gene_up <- intersect(which(gene$log2FoldChange >= log2(1.5)),which(gene$padj <= 0.05))
gene_down <- intersect(which(gene$log2FoldChange <= -log2(1.5)), which(gene$padj <= 0.05))
sig <- rep("no", times = nrow(gene))
sig[gene_up] <- "up"
sig[gene_down] <- "down"   
sig <- factor(sig,levels=c("up","down","no")) 
gene$type <- sig
gene$log10FDR <- -log(gene$padj, 15)
gene$log10FDR2 <- gene$log10FDR
gene$log10FDR2[gene$log10FDR2 > 30] <- 30
gene$log2FoldChange2 <- gene$log2FoldChange
gene$log2FoldChange2[gene$log2FoldChange2 > 15] <- 15
gene$log2FoldChange2[gene$log2FoldChange2 < -15] <- -15
#挑出上下调的去画，不然点太多，放在AI里面太卡
gene <- gene[which(gene$type=="up" | gene$type=="down"),]
#画图
library("ggplot2")
pdf("volcanoplot ehd3 vs WT padj0.05 log2(1.5).pdf",w = 6,h=4)
ggplot(gene,aes(x = log10FDR2 ,y=log2FoldChange2,col = gene$type))+
  geom_point(size = 0.5 ,alpha = 0.6) + 
  labs(x = "-log10(padj)", y =  "log2(fold change)") +theme_bw()+
  scale_colour_manual(values = c("up" = "red", "down" = "green"))+
  geom_hline(yintercept = c(-log2(1.5),log2(1.5)),linetype = "dashed",, colour = "grey") +
  geom_vline(xintercept = c(-log(0.05, 10)),linetype = "dashed", colour = "grey") +
  ylim(-15,15) + xlim(0,30) +
  theme_bw() + theme(panel.background = element_rect(colour = "black", size = 1, fill = "white"), panel.grid = element_blank())
dev.off()
length(gene_up);length(gene_down)
#1233 806

#volcanoplot sdg724 vs WT padj0.05 log2(1.5)
setwd(dir = "C:/Users/HP/Desktop")
gene <- read.delim("E:/复旦/实验室/课题/8_炳师兄水稻/new/RNAseq/14days/差异基因/合并计算/padj0.05 log2(1.5)/sdg724 vs WT padj0.05 log2(1.5).xls")
gene_up <- intersect(which(gene$log2FoldChange >= log2(1.5)),which(gene$padj <= 0.05))
gene_down <- intersect(which(gene$log2FoldChange <= -log2(1.5)), which(gene$padj <= 0.05))
sig <- rep("no", times = nrow(gene))
sig[gene_up] <- "up"
sig[gene_down] <- "down"   
sig <- factor(sig,levels=c("up","down","no")) 
gene$type <- sig
gene$log10FDR <- -log(gene$padj, 15)
gene$log10FDR2 <- gene$log10FDR
gene$log10FDR2[gene$log10FDR2 > 30] <- 30
gene$log2FoldChange2 <- gene$log2FoldChange
gene$log2FoldChange2[gene$log2FoldChange2 > 15] <- 15
gene$log2FoldChange2[gene$log2FoldChange2 < -15] <- -15
#挑出上下调的去画，不然点太多，放在AI里面太卡
gene <- gene[which(gene$type=="up" | gene$type=="down"),]
#画图
library("ggplot2")
pdf("volcanoplot sdg724 vs WT padj0.05 log2(1.5).pdf",w = 6,h=4)
ggplot(gene,aes(x = log10FDR2 ,y=log2FoldChange2,col = gene$type))+
  geom_point(size = 0.5 ,alpha = 0.6) + 
  labs(x = "-log10(padj)", y =  "log2(fold change)") +theme_bw()+
  scale_colour_manual(values = c("up" = "red", "down" = "green"))+
  geom_hline(yintercept = c(-log2(1.5),log2(1.5)),linetype = "dashed",, colour = "grey") +
  geom_vline(xintercept = c(-log(0.05, 10)),linetype = "dashed", colour = "grey") +
  ylim(-15,15) + xlim(0,30) +
  theme_bw() + theme(panel.background = element_rect(colour = "black", size = 1, fill = "white"), panel.grid = element_blank())
dev.off()
length(gene_up);length(gene_down)
#1694 1213


#volcanoplot ehd3sdg724 vs WT padj0.05 log2(1.5)
setwd(dir = "C:/Users/HP/Desktop")
gene <- read.delim("E:/复旦/实验室/课题/8_炳师兄水稻/new/RNAseq/14days/差异基因/合并计算/padj0.05 log2(1.5)/ehd3sdg724 vs WT padj0.05 log2(1.5).xls")
gene_up <- intersect(which(gene$log2FoldChange >= log2(1.5)),which(gene$padj <= 0.05))
gene_down <- intersect(which(gene$log2FoldChange <= -log2(1.5)), which(gene$padj <= 0.05))
sig <- rep("no", times = nrow(gene))
sig[gene_up] <- "up"
sig[gene_down] <- "down"   
sig <- factor(sig,levels=c("up","down","no")) 
gene$type <- sig
gene$log10FDR <- -log(gene$padj, 15)
gene$log10FDR2 <- gene$log10FDR
gene$log10FDR2[gene$log10FDR2 > 30] <- 30
gene$log2FoldChange2 <- gene$log2FoldChange
gene$log2FoldChange2[gene$log2FoldChange2 > 15] <- 15
gene$log2FoldChange2[gene$log2FoldChange2 < -15] <- -15
#挑出上下调的去画，不然点太多，放在AI里面太卡
gene <- gene[which(gene$type=="up" | gene$type=="down"),]
#画图
library("ggplot2")
pdf("volcanoplot ehd3sdg724 vs WT padj0.05 log2(1.5).pdf",w = 6,h=4)
ggplot(gene,aes(x = log10FDR2 ,y=log2FoldChange2,col = gene$type))+
  geom_point(size = 0.5 ,alpha = 0.6) + 
  labs(x = "-log10(padj)", y =  "log2(fold change)") +theme_bw()+
  scale_colour_manual(values = c("up" = "red", "down" = "green"))+
  geom_hline(yintercept = c(-log2(1.5),log2(1.5)),linetype = "dashed",, colour = "grey") +
  geom_vline(xintercept = c(-log(0.05, 10)),linetype = "dashed", colour = "grey") +
  ylim(-15,15) + xlim(0,30) +
  theme_bw() + theme(panel.background = element_rect(colour = "black", size = 1, fill = "white"), panel.grid = element_blank())
dev.off()
length(gene_up);length(gene_down)
#951  704














#GO analysis
suppressMessages(library(clusterProfiler))
suppressMessages(library(topGO))
suppressMessages(library(org.At.tair.db))
genelist <- read.delim("fas2_up.txt",header=F)
eGOBP <- enrichGO(gene = genelist[,1], OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1, qvalueCutoff = 1, readable = F)
eGOBP1 <- enrichGO(gene = genelist[,1], OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = F)
write.table(data.frame(eGOBP),"fas2_up GO Enrichment_BP.xls",sep="\t",row.names=F)
barplot(eGOBP1, showCategory=20,title="Enrichment_BP")
#dotplot(eGOBP1, showCategory=20, title="Enrichment_BP")

#其实它R包自带的画图很弱，推荐使用它生成的excel自己去网站画图，或者有问题我们再交流
#排序、筛选、根据p还是q还是别的，有待自行探索
#还有很多在线工具可以做，比如水稻的GO是用一个好像叫carmo的做的，maybe，拟南芥还可以用David GO做，也是个在线网站

