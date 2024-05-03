#DEGs
setwd("C:/Users/94526/Desktop/")
#If the number of replicates is not 3, change the following section on Replicates to a different number.
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
#############################################
sp1  <- read.delim("WT_rep1.rawcount")
sp2  <- read.delim("WT_rep2.rawcount")
sp3  <- read.delim("WT_rep3.rawcount")
sp4  <- read.delim("mutant_rep1.rawcount")
sp5  <- read.delim("mutant_rep2.rawcount")
sp6  <- read.delim("mutant_rep3.rawcount")
DEGfind("WT","mutant","1.5","2","0.05","0.01")
