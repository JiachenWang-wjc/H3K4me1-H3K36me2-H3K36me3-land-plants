setwd("/media/wangjiachen/disk2/workflow/ehd3/00_all_cellection_process/heatmap_fuke_xiewenhao")
allm <- read.delim("all_modification_n55801.mat", header = F, skip=1)
row.names(allm) <- allm[,4]
nrow(allm)
ncol(allm)
H3K4me1 <- allm[, 7:506]
H3K4me2 <- allm[, 507:1006]
H3K4me3 <- allm[, 1007:1506]
H3K36me1 <- allm[, 1507:2006]
H3K36me2 <- allm[, 2007:2506]
H3K36me3 <- allm[, 2507:3006]
H3K9me1 <- allm[, 3007:3506]
H3K9me2 <- allm[, 3507:4006]
H3K9me3 <- allm[, 4007:4506]
H3K27me3 <- allm[, 4507:5006]
H3K4ac <- allm[, 5007:5506]
H3K9ac <- allm[, 5507:6006]
H3K27ac <- allm[, 6007:6506]
H3K56ac <- allm[, 6507:7006]
H4K12ac <- allm[, 7007:7506]

H3K4me1_GB <- allm[, 107:406]
H3K4me2_GB <- allm[, 607:906]
H3K4me3_GB <- allm[, 1107:1406]
H3K36me1_GB <- allm[, 1607:1906]
H3K36me2_GB <- allm[, 2107:2406]
H3K36me3_GB <- allm[, 2607:2906]
H3K9me1_GB <- allm[, 3107:3406]
H3K9me2_GB <- allm[, 3607:3906]
H3K9me3_GB <- allm[, 4107:4406]
H3K27me3_GB <- allm[, 4607:4906]
H3K4ac_GB <- allm[, 5107:5406]
H3K9ac_GB <- allm[, 5607:5906]
H3K27ac_GB <- allm[, 6107:6406]
H3K56ac_GB <- allm[, 6607:6906]
H4K12ac_GB <- allm[, 7107:7406]





rowMeanDataList <- list(rowMeans(H3K4me1_GB), rowMeans(H3K4me2_GB), rowMeans(H3K4me3_GB), rowMeans(H3K36me1_GB), rowMeans(H3K36me2_GB), rowMeans(H3K36me3_GB), rowMeans(H3K9me1_GB), rowMeans(H3K9me2_GB), rowMeans(H3K9me3_GB), rowMeans(H3K27me3_GB), rowMeans(H3K4ac_GB), rowMeans(H3K9ac_GB), rowMeans(H3K27ac_GB), rowMeans(H3K56ac_GB), rowMeans(H4K12ac_GB))
write.table(rowMeanDataList, "GB_modification_level.txt", col.names = F, row.names = T,sep = "\t",quote = F)
write.table(rowMeanDataList, "GB_modification_level.xls", col.names = F, row.names = T,sep = "\t",quote = F)

#用IGV检查过了，没什么问题，得到一个txt和一个xls



GB_modi_level <- read.delim("GB_modification_level.txt", header = F)
row.names(GB_modi_level) <- GB_modi_level[,1]
GB_modi_level_new <- GB_modi_level[,2:16]
nrow(GB_modi_level_new)
ncol(GB_modi_level_new)
#55801
#15

#write.table(GB_modi_level_new, "test.xls", col.names = F, row.names = T,sep = "\t",quote = F)
#没问题，这下把第一列基因名变成行名了


#每一列是一种修饰。每一行都是一个基因











#数据框有55801行，15列，求每一列的最小值，最大值，最大值和最小值的差值，均值,标准差，然后进行标准化
#找出每一列的最小值
min <- apply(GB_modi_level_new, MARGIN = 2, FUN = min)
#找出每一列的最大值
max <- apply(GB_modi_level_new, MARGIN = 2, FUN = max) 
#求最大值和最小值的差值
diff <- max - min
#求每一列的均值
mean <- apply(GB_modi_level_new, MARGIN = 2, FUN = mean) 
#求每一列的标准差
sd <- apply(GB_modi_level_new, MARGIN = 2, FUN = sd) 




#每个值减去减去均值，然后再除以那一行的标准差，就是z-normalizing后的数值
GB_modi_level_new$H3K4me1 <- (GB_modi_level_new[,1]-mean[[1]]) / sd[[1]]
GB_modi_level_new$H3K4me2 <- (GB_modi_level_new[,2]-mean[[2]]) / sd[[2]]
GB_modi_level_new$H3K4me3 <- (GB_modi_level_new[,3]-mean[[3]]) / sd[[3]]
GB_modi_level_new$H3K36me1 <- (GB_modi_level_new[,4]-mean[[4]]) / sd[[4]]
GB_modi_level_new$H3K36me2 <- (GB_modi_level_new[,5]-mean[[5]]) / sd[[5]]
GB_modi_level_new$H3K36me3 <- (GB_modi_level_new[,6]-mean[[6]]) / sd[[6]]
GB_modi_level_new$H3K9me1 <- (GB_modi_level_new[,7]-mean[[7]]) / sd[[7]]
GB_modi_level_new$H3K9me2 <- (GB_modi_level_new[,8]-mean[[8]]) / sd[[8]]
GB_modi_level_new$H3K9me3 <- (GB_modi_level_new[,9]-mean[[9]]) / sd[[9]]
GB_modi_level_new$H3K27me3 <- (GB_modi_level_new[,10]-mean[[10]]) / sd[[10]]
GB_modi_level_new$H3K4ac <- (GB_modi_level_new[,11]-mean[[11]]) / sd[[11]]
GB_modi_level_new$H3K9ac <- (GB_modi_level_new[,12]-mean[[12]]) / sd[[12]]
GB_modi_level_new$H3K27ac <- (GB_modi_level_new[,13]-mean[[13]]) / sd[[13]]
GB_modi_level_new$H3K56ac <- (GB_modi_level_new[,14]-mean[[14]]) / sd[[14]]
GB_modi_level_new$H4K12ac <- (GB_modi_level_new[,15]-mean[[15]]) / sd[[15]]
#保存数据
write.table(GB_modi_level_new, file = "Allgenes_modification_z-normalized.txt", col.names = F, row.names = T,sep = "\t",quote = F)
write.table(GB_modi_level_new, file = "Allgenes_modification_z-normalized.xls", col.names = F, row.names = T,sep = "\t",quote = F)
#通过IGV看也是没什么问题的





#下面用我自己的电脑画图，服务器上面貌似装不了那些R包




setwd("F:/000实验室/111在进行的课题/ehd3/00_投第一稿前的大整合重新画图/各种尝试/复刻谢文浩不同修饰相关性热图")
All_zdata <- read.delim("Allgenes_modification_z-normalized.xls", header = F)
nrow(All_zdata)
#55801
ncol(All_zdata)
#31
row.names(All_zdata) <- All_zdata[,1]
All_zdata_new <- All_zdata[,17:31]
nrow(All_zdata_new)
ncol(All_zdata_new)
#15
library(ComplexHeatmap)  #提供Heatmap函数
library(circlize)  #提供colorRamp2函数
library(RColorBrewer) #提供brewer.pal函数
heatdata <- All_zdata_new
k=kmeans(heatdata,50)
heatdata$k=k$cluster
#根据K值对heatdatas数据重新排序
htframe_sort=heatdata[order(heatdata$k),]
#将数据框转化成矩阵，因为之后heatmap()函数需要输入的是矩阵
htmatrix_sort <- as.matrix(htframe_sort[,1:14])
#heatmap绘图,show_row_names：是否显示行名称。默认值为TRUE
#show_column_names：是否显示列名称。默认值为TRUE
#cluster_rows=FALSE,如果为TRUE，则在行上创建聚类的簇；
#cluster_columns=FALSE,如果为TRUE，则在列上创建聚类的簇。
write.table(htframe_sort,"Allgenes_modification_z-normalized_heatmap1.xls",sep="\t")






All_zdata <- read.delim("Allgenes_modification_z-normalized.xls", header = F)
nrow(All_zdata)
#55801
ncol(All_zdata)
#31
row.names(All_zdata) <- All_zdata[,1]
All_zdata_new <- All_zdata[,17:31]
nrow(All_zdata_new)
ncol(All_zdata_new)






##计算pearson相关系数不会报错，计算spearman会报错： Cannot compute exact p-value with ties，但可以输出结果

#计算spearman相关系数 H3K4me1和其他的marker,并将结果放到同一个向量
H3K4me1_H3K4me1  <- cor.test(All_zdata_new[,1], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K4me1_H3K4me2  <- cor.test(All_zdata_new[,1], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K4me1_H3K4me3  <- cor.test(All_zdata_new[,1], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K4me1_H3K36me1 <- cor.test(All_zdata_new[,1], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K4me1_H3K36me2 <- cor.test(All_zdata_new[,1], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K4me1_H3K36me3 <- cor.test(All_zdata_new[,1], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9me1   <- cor.test(All_zdata_new[,1], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9me2   <- cor.test(All_zdata_new[,1], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9me3   <- cor.test(All_zdata_new[,1], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K4me1_H3K27me3   <- cor.test(All_zdata_new[,1], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K4me1_H3K4ac   <- cor.test(All_zdata_new[,1], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9ac   <- cor.test(All_zdata_new[,1], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K4me1_H3K27ac  <- cor.test(All_zdata_new[,1], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K4me1_H3K56ac  <- cor.test(All_zdata_new[,1], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K4me1_H4K12ac   <- cor.test(All_zdata_new[,1], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x1 <- c(H3K4me1_H3K4me1$estimate, H3K4me1_H3K4me2$estimate, H3K4me1_H3K4me3$estimate, H3K4me1_H3K36me1$estimate,H3K4me1_H3K36me2$estimate, H3K4me1_H3K36me3$estimate, 
        H3K4me1_H3K9me1$estimate, H3K4me1_H3K9me2$estimate, H3K4me1_H3K9me3$estimate,H3K4me1_H3K27me3$estimate,
        H3K4me1_H3K4ac$estimate, H3K4me1_H3K9ac$estimate, H3K4me1_H3K27ac$estimate, H3K4me1_H3K56ac$estimate, H3K4me1_H4K12ac$estimate)

#结果不是很理想，还是按照只有H3K4me1的基因作为背景画图试一下看能不能和谢文吻合





All_zdata <- read.delim("Allgenes_modification_z-normalized.xls", header = F)
nrow(All_zdata)
#55801
ncol(All_zdata)
#31
row.names(All_zdata) <- All_zdata[,1]
All_zdata_new <- All_zdata[,17:31]
nrow(All_zdata_new)
ncol(All_zdata_new)
#55801
#15






target <- read.delim("F:/000实验室/111在进行的课题/ehd3/00_投第一稿前的大整合重新画图/各种尝试/倒着看自身定位/H3K4me1的各种lists/Batch2_spikein_WT_H3K4me1_q005_broad_peakRelatedGenes_annoOvgene_genelist.xls")
row.names(target) <- target[,1]
All_zdata_new <- subset(All_zdata_new, is.element(row.names(All_zdata_new), row.names(target)) == T)
#19323
#15

#计算spearman相关系数 H3K4me1和其他的marker,并将结果放到同一个向量
H3K4me1_H3K4me1  <- cor.test(All_zdata_new[,1], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K4me1_H3K4me2  <- cor.test(All_zdata_new[,1], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K4me1_H3K4me3  <- cor.test(All_zdata_new[,1], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K4me1_H3K36me1 <- cor.test(All_zdata_new[,1], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K4me1_H3K36me2 <- cor.test(All_zdata_new[,1], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K4me1_H3K36me3 <- cor.test(All_zdata_new[,1], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9me1   <- cor.test(All_zdata_new[,1], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9me2   <- cor.test(All_zdata_new[,1], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9me3   <- cor.test(All_zdata_new[,1], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K4me1_H3K27me3   <- cor.test(All_zdata_new[,1], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K4me1_H3K4ac   <- cor.test(All_zdata_new[,1], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K4me1_H3K9ac   <- cor.test(All_zdata_new[,1], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K4me1_H3K27ac  <- cor.test(All_zdata_new[,1], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K4me1_H3K56ac  <- cor.test(All_zdata_new[,1], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K4me1_H4K12ac   <- cor.test(All_zdata_new[,1], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x1 <- c(H3K4me1_H3K4me1$estimate, H3K4me1_H3K4me2$estimate, H3K4me1_H3K4me3$estimate, H3K4me1_H3K36me1$estimate,H3K4me1_H3K36me2$estimate, H3K4me1_H3K36me3$estimate, 
        H3K4me1_H3K9me1$estimate, H3K4me1_H3K9me2$estimate, H3K4me1_H3K9me3$estimate,H3K4me1_H3K27me3$estimate,
        H3K4me1_H3K4ac$estimate, H3K4me1_H3K9ac$estimate, H3K4me1_H3K27ac$estimate, H3K4me1_H3K56ac$estimate, H3K4me1_H4K12ac$estimate)

#结果不是很理想，还是按照只有H3K4me1的基因作为背景画图试一下看能不能和谢文吻合
#霍好家伙，一下子变得吻合了，艹，看来全基因还是不行啊



H3K4me2_H3K4me1  <- cor.test(All_zdata_new[,2], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K4me2_H3K4me2  <- cor.test(All_zdata_new[,2], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K4me2_H3K4me3  <- cor.test(All_zdata_new[,2], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K4me2_H3K36me1 <- cor.test(All_zdata_new[,2], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K4me2_H3K36me2 <- cor.test(All_zdata_new[,2], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K4me2_H3K36me3 <- cor.test(All_zdata_new[,2], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K4me2_H3K9me1   <- cor.test(All_zdata_new[,2], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K4me2_H3K9me2   <- cor.test(All_zdata_new[,2], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K4me2_H3K9me3   <- cor.test(All_zdata_new[,2], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K4me2_H3K27me3   <- cor.test(All_zdata_new[,2], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K4me2_H3K4ac   <- cor.test(All_zdata_new[,2], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K4me2_H3K9ac   <- cor.test(All_zdata_new[,2], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K4me2_H3K27ac  <- cor.test(All_zdata_new[,2], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K4me2_H3K56ac  <- cor.test(All_zdata_new[,2], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K4me2_H4K12ac   <- cor.test(All_zdata_new[,2], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x2 <- c(H3K4me2_H3K4me1$estimate, H3K4me2_H3K4me2$estimate, H3K4me2_H3K4me3$estimate, H3K4me2_H3K36me1$estimate,H3K4me2_H3K36me2$estimate, H3K4me2_H3K36me3$estimate, 
        H3K4me2_H3K9me1$estimate, H3K4me2_H3K9me2$estimate, H3K4me2_H3K9me3$estimate,H3K4me2_H3K27me3$estimate,
        H3K4me2_H3K4ac$estimate, H3K4me2_H3K9ac$estimate, H3K4me2_H3K27ac$estimate, H3K4me2_H3K56ac$estimate, H3K4me2_H4K12ac$estimate)





H3K4me3_H3K4me1  <- cor.test(All_zdata_new[,3], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K4me3_H3K4me2  <- cor.test(All_zdata_new[,3], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K4me3_H3K4me3  <- cor.test(All_zdata_new[,3], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K4me3_H3K36me1 <- cor.test(All_zdata_new[,3], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K4me3_H3K36me2 <- cor.test(All_zdata_new[,3], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K4me3_H3K36me3 <- cor.test(All_zdata_new[,3], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K4me3_H3K9me1   <- cor.test(All_zdata_new[,3], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K4me3_H3K9me2   <- cor.test(All_zdata_new[,3], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K4me3_H3K9me3   <- cor.test(All_zdata_new[,3], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K4me3_H3K27me3   <- cor.test(All_zdata_new[,3], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K4me3_H3K4ac   <- cor.test(All_zdata_new[,3], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K4me3_H3K9ac   <- cor.test(All_zdata_new[,3], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K4me3_H3K27ac  <- cor.test(All_zdata_new[,3], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K4me3_H3K56ac  <- cor.test(All_zdata_new[,3], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K4me3_H4K12ac   <- cor.test(All_zdata_new[,3], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x3 <- c(H3K4me3_H3K4me1$estimate, H3K4me3_H3K4me2$estimate, H3K4me3_H3K4me3$estimate, H3K4me3_H3K36me1$estimate,H3K4me3_H3K36me2$estimate, H3K4me3_H3K36me3$estimate, 
        H3K4me3_H3K9me1$estimate, H3K4me3_H3K9me2$estimate, H3K4me3_H3K9me3$estimate,H3K4me3_H3K27me3$estimate,
        H3K4me3_H3K4ac$estimate, H3K4me3_H3K9ac$estimate, H3K4me3_H3K27ac$estimate, H3K4me3_H3K56ac$estimate, H3K4me3_H4K12ac$estimate)





H3K36me1_H3K4me1  <- cor.test(All_zdata_new[,4], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K36me1_H3K4me2  <- cor.test(All_zdata_new[,4], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K36me1_H3K4me3  <- cor.test(All_zdata_new[,4], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K36me1_H3K36me1 <- cor.test(All_zdata_new[,4], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K36me1_H3K36me2 <- cor.test(All_zdata_new[,4], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K36me1_H3K36me3 <- cor.test(All_zdata_new[,4], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K36me1_H3K9me1   <- cor.test(All_zdata_new[,4], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K36me1_H3K9me2   <- cor.test(All_zdata_new[,4], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K36me1_H3K9me3   <- cor.test(All_zdata_new[,4], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K36me1_H3K27me3   <- cor.test(All_zdata_new[,4], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K36me1_H3K4ac   <- cor.test(All_zdata_new[,4], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K36me1_H3K9ac   <- cor.test(All_zdata_new[,4], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K36me1_H3K27ac  <- cor.test(All_zdata_new[,4], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K36me1_H3K56ac  <- cor.test(All_zdata_new[,4], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K36me1_H4K12ac   <- cor.test(All_zdata_new[,4], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x4 <- c(H3K36me1_H3K4me1$estimate, H3K36me1_H3K4me2$estimate, H3K36me1_H3K4me3$estimate, H3K36me1_H3K36me1$estimate,H3K36me1_H3K36me2$estimate, H3K36me1_H3K36me3$estimate, 
        H3K36me1_H3K9me1$estimate, H3K36me1_H3K9me2$estimate, H3K36me1_H3K9me3$estimate,H3K36me1_H3K27me3$estimate,
        H3K36me1_H3K4ac$estimate, H3K36me1_H3K9ac$estimate, H3K36me1_H3K27ac$estimate, H3K36me1_H3K56ac$estimate, H3K36me1_H4K12ac$estimate)



H3K36me2_H3K4me1  <- cor.test(All_zdata_new[,5], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K36me2_H3K4me2  <- cor.test(All_zdata_new[,5], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K36me2_H3K4me3  <- cor.test(All_zdata_new[,5], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K36me2_H3K36me1 <- cor.test(All_zdata_new[,5], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K36me2_H3K36me2 <- cor.test(All_zdata_new[,5], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K36me2_H3K36me3 <- cor.test(All_zdata_new[,5], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K36me2_H3K9me1   <- cor.test(All_zdata_new[,5], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K36me2_H3K9me2   <- cor.test(All_zdata_new[,5], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K36me2_H3K9me3   <- cor.test(All_zdata_new[,5], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K36me2_H3K27me3   <- cor.test(All_zdata_new[,5], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K36me2_H3K4ac   <- cor.test(All_zdata_new[,5], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K36me2_H3K9ac   <- cor.test(All_zdata_new[,5], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K36me2_H3K27ac  <- cor.test(All_zdata_new[,5], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K36me2_H3K56ac  <- cor.test(All_zdata_new[,5], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K36me2_H4K12ac   <- cor.test(All_zdata_new[,5], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x5 <- c(H3K36me2_H3K4me1$estimate, H3K36me2_H3K4me2$estimate, H3K36me2_H3K4me3$estimate, H3K36me2_H3K36me1$estimate,H3K36me2_H3K36me2$estimate, H3K36me2_H3K36me3$estimate, 
        H3K36me2_H3K9me1$estimate, H3K36me2_H3K9me2$estimate, H3K36me2_H3K9me3$estimate,H3K36me2_H3K27me3$estimate,
        H3K36me2_H3K4ac$estimate, H3K36me2_H3K9ac$estimate, H3K36me2_H3K27ac$estimate, H3K36me2_H3K56ac$estimate, H3K36me2_H4K12ac$estimate)


        
H3K36me3_H3K4me1  <- cor.test(All_zdata_new[,6], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K36me3_H3K4me2  <- cor.test(All_zdata_new[,6], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K36me3_H3K4me3  <- cor.test(All_zdata_new[,6], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K36me3_H3K36me1 <- cor.test(All_zdata_new[,6], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K36me3_H3K36me2 <- cor.test(All_zdata_new[,6], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K36me3_H3K36me3 <- cor.test(All_zdata_new[,6], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K36me3_H3K9me1   <- cor.test(All_zdata_new[,6], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K36me3_H3K9me2   <- cor.test(All_zdata_new[,6], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K36me3_H3K9me3   <- cor.test(All_zdata_new[,6], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K36me3_H3K27me3   <- cor.test(All_zdata_new[,6], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K36me3_H3K4ac   <- cor.test(All_zdata_new[,6], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K36me3_H3K9ac   <- cor.test(All_zdata_new[,6], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K36me3_H3K27ac  <- cor.test(All_zdata_new[,6], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K36me3_H3K56ac  <- cor.test(All_zdata_new[,6], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K36me3_H4K12ac   <- cor.test(All_zdata_new[,6], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x6 <- c(H3K36me3_H3K4me1$estimate, H3K36me3_H3K4me2$estimate, H3K36me3_H3K4me3$estimate, H3K36me3_H3K36me1$estimate,H3K36me3_H3K36me2$estimate, H3K36me3_H3K36me3$estimate, 
        H3K36me3_H3K9me1$estimate, H3K36me3_H3K9me2$estimate, H3K36me3_H3K9me3$estimate,H3K36me3_H3K27me3$estimate,
        H3K36me3_H3K4ac$estimate, H3K36me3_H3K9ac$estimate, H3K36me3_H3K27ac$estimate, H3K36me3_H3K56ac$estimate, H3K36me3_H4K12ac$estimate)





H3K9me1_H3K4me1  <- cor.test(All_zdata_new[,7], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K9me1_H3K4me2  <- cor.test(All_zdata_new[,7], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K9me1_H3K4me3  <- cor.test(All_zdata_new[,7], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K9me1_H3K36me1 <- cor.test(All_zdata_new[,7], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K9me1_H3K36me2 <- cor.test(All_zdata_new[,7], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K9me1_H3K36me3 <- cor.test(All_zdata_new[,7], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K9me1_H3K9me1   <- cor.test(All_zdata_new[,7], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K9me1_H3K9me2   <- cor.test(All_zdata_new[,7], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K9me1_H3K9me3   <- cor.test(All_zdata_new[,7], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K9me1_H3K27me3   <- cor.test(All_zdata_new[,7], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K9me1_H3K4ac   <- cor.test(All_zdata_new[,7], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K9me1_H3K9ac   <- cor.test(All_zdata_new[,7], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K9me1_H3K27ac  <- cor.test(All_zdata_new[,7], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K9me1_H3K56ac  <- cor.test(All_zdata_new[,7], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K9me1_H4K12ac   <- cor.test(All_zdata_new[,7], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x7 <- c(H3K9me1_H3K4me1$estimate, H3K9me1_H3K4me2$estimate, H3K9me1_H3K4me3$estimate, H3K9me1_H3K36me1$estimate,H3K9me1_H3K36me2$estimate, H3K9me1_H3K36me3$estimate, 
        H3K9me1_H3K9me1$estimate, H3K9me1_H3K9me2$estimate, H3K9me1_H3K9me3$estimate,H3K9me1_H3K27me3$estimate,
        H3K9me1_H3K4ac$estimate, H3K9me1_H3K9ac$estimate, H3K9me1_H3K27ac$estimate, H3K9me1_H3K56ac$estimate, H3K9me1_H4K12ac$estimate)







H3K9me2_H3K4me1  <- cor.test(All_zdata_new[,8], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K9me2_H3K4me2  <- cor.test(All_zdata_new[,8], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K9me2_H3K4me3  <- cor.test(All_zdata_new[,8], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K9me2_H3K36me1 <- cor.test(All_zdata_new[,8], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K9me2_H3K36me2 <- cor.test(All_zdata_new[,8], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K9me2_H3K36me3 <- cor.test(All_zdata_new[,8], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K9me2_H3K9me1   <- cor.test(All_zdata_new[,8], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K9me2_H3K9me2   <- cor.test(All_zdata_new[,8], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K9me2_H3K9me3   <- cor.test(All_zdata_new[,8], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K9me2_H3K27me3   <- cor.test(All_zdata_new[,8], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K9me2_H3K4ac   <- cor.test(All_zdata_new[,8], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K9me2_H3K9ac   <- cor.test(All_zdata_new[,8], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K9me2_H3K27ac  <- cor.test(All_zdata_new[,8], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K9me2_H3K56ac  <- cor.test(All_zdata_new[,8], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K9me2_H4K12ac   <- cor.test(All_zdata_new[,8], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x8 <- c(H3K9me2_H3K4me1$estimate, H3K9me2_H3K4me2$estimate, H3K9me2_H3K4me3$estimate, H3K9me2_H3K36me1$estimate,H3K9me2_H3K36me2$estimate, H3K9me2_H3K36me3$estimate, 
        H3K9me2_H3K9me1$estimate, H3K9me2_H3K9me2$estimate, H3K9me2_H3K9me3$estimate,H3K9me2_H3K27me3$estimate,
        H3K9me2_H3K4ac$estimate, H3K9me2_H3K9ac$estimate, H3K9me2_H3K27ac$estimate, H3K9me2_H3K56ac$estimate, H3K9me2_H4K12ac$estimate)






H3K9me3_H3K4me1  <- cor.test(All_zdata_new[,9], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K9me3_H3K4me2  <- cor.test(All_zdata_new[,9], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K9me3_H3K4me3  <- cor.test(All_zdata_new[,9], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K9me3_H3K36me1 <- cor.test(All_zdata_new[,9], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K9me3_H3K36me2 <- cor.test(All_zdata_new[,9], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K9me3_H3K36me3 <- cor.test(All_zdata_new[,9], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K9me3_H3K9me1   <- cor.test(All_zdata_new[,9], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K9me3_H3K9me2   <- cor.test(All_zdata_new[,9], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K9me3_H3K9me3   <- cor.test(All_zdata_new[,9], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K9me3_H3K27me3   <- cor.test(All_zdata_new[,9], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K9me3_H3K4ac   <- cor.test(All_zdata_new[,9], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K9me3_H3K9ac   <- cor.test(All_zdata_new[,9], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K9me3_H3K27ac  <- cor.test(All_zdata_new[,9], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K9me3_H3K56ac  <- cor.test(All_zdata_new[,9], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K9me3_H4K12ac   <- cor.test(All_zdata_new[,9], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x9 <- c(H3K9me3_H3K4me1$estimate, H3K9me3_H3K4me2$estimate, H3K9me3_H3K4me3$estimate, H3K9me3_H3K36me1$estimate,H3K9me3_H3K36me2$estimate, H3K9me3_H3K36me3$estimate, 
        H3K9me3_H3K9me1$estimate, H3K9me3_H3K9me2$estimate, H3K9me3_H3K9me3$estimate,H3K9me3_H3K27me3$estimate,
        H3K9me3_H3K4ac$estimate, H3K9me3_H3K9ac$estimate, H3K9me3_H3K27ac$estimate, H3K9me3_H3K56ac$estimate, H3K9me3_H4K12ac$estimate)





H3K27me3_H3K4me1  <- cor.test(All_zdata_new[,10], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K27me3_H3K4me2  <- cor.test(All_zdata_new[,10], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K27me3_H3K4me3  <- cor.test(All_zdata_new[,10], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K27me3_H3K36me1 <- cor.test(All_zdata_new[,10], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K27me3_H3K36me2 <- cor.test(All_zdata_new[,10], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K27me3_H3K36me3 <- cor.test(All_zdata_new[,10], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K27me3_H3K9me1   <- cor.test(All_zdata_new[,10], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K27me3_H3K9me2   <- cor.test(All_zdata_new[,10], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K27me3_H3K9me3   <- cor.test(All_zdata_new[,10], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K27me3_H3K27me3   <- cor.test(All_zdata_new[,10], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K27me3_H3K4ac   <- cor.test(All_zdata_new[,10], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K27me3_H3K9ac   <- cor.test(All_zdata_new[,10], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K27me3_H3K27ac  <- cor.test(All_zdata_new[,10], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K27me3_H3K56ac  <- cor.test(All_zdata_new[,10], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K27me3_H4K12ac   <- cor.test(All_zdata_new[,10], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x10 <- c(H3K27me3_H3K4me1$estimate, H3K27me3_H3K4me2$estimate, H3K27me3_H3K4me3$estimate, H3K27me3_H3K36me1$estimate,H3K27me3_H3K36me2$estimate, H3K27me3_H3K36me3$estimate, 
        H3K27me3_H3K9me1$estimate, H3K27me3_H3K9me2$estimate, H3K27me3_H3K9me3$estimate,H3K27me3_H3K27me3$estimate,
        H3K27me3_H3K4ac$estimate, H3K27me3_H3K9ac$estimate, H3K27me3_H3K27ac$estimate, H3K27me3_H3K56ac$estimate, H3K27me3_H4K12ac$estimate)






H3K4ac_H3K4me1  <- cor.test(All_zdata_new[,11], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K4ac_H3K4me2  <- cor.test(All_zdata_new[,11], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K4ac_H3K4me3  <- cor.test(All_zdata_new[,11], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K4ac_H3K36me1 <- cor.test(All_zdata_new[,11], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K4ac_H3K36me2 <- cor.test(All_zdata_new[,11], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K4ac_H3K36me3 <- cor.test(All_zdata_new[,11], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K4ac_H3K9me1   <- cor.test(All_zdata_new[,11], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K4ac_H3K9me2   <- cor.test(All_zdata_new[,11], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K4ac_H3K9me3   <- cor.test(All_zdata_new[,11], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K4ac_H3K27me3   <- cor.test(All_zdata_new[,11], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K4ac_H3K4ac   <- cor.test(All_zdata_new[,11], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K4ac_H3K9ac   <- cor.test(All_zdata_new[,11], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K4ac_H3K27ac  <- cor.test(All_zdata_new[,11], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K4ac_H3K56ac  <- cor.test(All_zdata_new[,11], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K4ac_H4K12ac   <- cor.test(All_zdata_new[,11], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x11 <- c(H3K4ac_H3K4me1$estimate, H3K4ac_H3K4me2$estimate, H3K4ac_H3K4me3$estimate, H3K4ac_H3K36me1$estimate,H3K4ac_H3K36me2$estimate, H3K4ac_H3K36me3$estimate, 
        H3K4ac_H3K9me1$estimate, H3K4ac_H3K9me2$estimate, H3K4ac_H3K9me3$estimate,H3K4ac_H3K27me3$estimate,
        H3K4ac_H3K4ac$estimate, H3K4ac_H3K9ac$estimate, H3K4ac_H3K27ac$estimate, H3K4ac_H3K56ac$estimate, H3K4ac_H4K12ac$estimate)









H3K9ac_H3K4me1  <- cor.test(All_zdata_new[,12], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K9ac_H3K4me2  <- cor.test(All_zdata_new[,12], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K9ac_H3K4me3  <- cor.test(All_zdata_new[,12], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K9ac_H3K36me1 <- cor.test(All_zdata_new[,12], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K9ac_H3K36me2 <- cor.test(All_zdata_new[,12], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K9ac_H3K36me3 <- cor.test(All_zdata_new[,12], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K9ac_H3K9me1   <- cor.test(All_zdata_new[,12], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K9ac_H3K9me2   <- cor.test(All_zdata_new[,12], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K9ac_H3K9me3   <- cor.test(All_zdata_new[,12], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K9ac_H3K27me3   <- cor.test(All_zdata_new[,12], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K9ac_H3K4ac   <- cor.test(All_zdata_new[,12], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K9ac_H3K9ac   <- cor.test(All_zdata_new[,12], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K9ac_H3K27ac  <- cor.test(All_zdata_new[,12], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K9ac_H3K56ac  <- cor.test(All_zdata_new[,12], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K9ac_H4K12ac   <- cor.test(All_zdata_new[,12], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x12 <- c(H3K9ac_H3K4me1$estimate, H3K9ac_H3K4me2$estimate, H3K9ac_H3K4me3$estimate, H3K9ac_H3K36me1$estimate,H3K9ac_H3K36me2$estimate, H3K9ac_H3K36me3$estimate, 
        H3K9ac_H3K9me1$estimate, H3K9ac_H3K9me2$estimate, H3K9ac_H3K9me3$estimate,H3K9ac_H3K27me3$estimate,
        H3K9ac_H3K4ac$estimate, H3K9ac_H3K9ac$estimate, H3K9ac_H3K27ac$estimate, H3K9ac_H3K56ac$estimate, H3K9ac_H4K12ac$estimate)





H3K27ac_H3K4me1  <- cor.test(All_zdata_new[,13], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K27ac_H3K4me2  <- cor.test(All_zdata_new[,13], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K27ac_H3K4me3  <- cor.test(All_zdata_new[,13], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K27ac_H3K36me1 <- cor.test(All_zdata_new[,13], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K27ac_H3K36me2 <- cor.test(All_zdata_new[,13], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K27ac_H3K36me3 <- cor.test(All_zdata_new[,13], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K27ac_H3K9me1   <- cor.test(All_zdata_new[,13], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K27ac_H3K9me2   <- cor.test(All_zdata_new[,13], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K27ac_H3K9me3   <- cor.test(All_zdata_new[,13], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K27ac_H3K27me3   <- cor.test(All_zdata_new[,13], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K27ac_H3K4ac   <- cor.test(All_zdata_new[,13], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K27ac_H3K9ac   <- cor.test(All_zdata_new[,13], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K27ac_H3K27ac  <- cor.test(All_zdata_new[,13], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K27ac_H3K56ac  <- cor.test(All_zdata_new[,13], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K27ac_H4K12ac   <- cor.test(All_zdata_new[,13], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x13 <- c(H3K27ac_H3K4me1$estimate, H3K27ac_H3K4me2$estimate, H3K27ac_H3K4me3$estimate, H3K27ac_H3K36me1$estimate,H3K27ac_H3K36me2$estimate, H3K27ac_H3K36me3$estimate, 
        H3K27ac_H3K9me1$estimate, H3K27ac_H3K9me2$estimate, H3K27ac_H3K9me3$estimate,H3K27ac_H3K27me3$estimate,
        H3K27ac_H3K4ac$estimate, H3K27ac_H3K9ac$estimate, H3K27ac_H3K27ac$estimate, H3K27ac_H3K56ac$estimate, H3K27ac_H4K12ac$estimate)





H3K56ac_H3K4me1  <- cor.test(All_zdata_new[,14], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H3K56ac_H3K4me2  <- cor.test(All_zdata_new[,14], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H3K56ac_H3K4me3  <- cor.test(All_zdata_new[,14], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H3K56ac_H3K36me1 <- cor.test(All_zdata_new[,14], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H3K56ac_H3K36me2 <- cor.test(All_zdata_new[,14], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H3K56ac_H3K36me3 <- cor.test(All_zdata_new[,14], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H3K56ac_H3K9me1   <- cor.test(All_zdata_new[,14], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H3K56ac_H3K9me2   <- cor.test(All_zdata_new[,14], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H3K56ac_H3K9me3   <- cor.test(All_zdata_new[,14], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H3K56ac_H3K27me3   <- cor.test(All_zdata_new[,14], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H3K56ac_H3K4ac   <- cor.test(All_zdata_new[,14], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H3K56ac_H3K9ac   <- cor.test(All_zdata_new[,14], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H3K56ac_H3K27ac  <- cor.test(All_zdata_new[,14], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H3K56ac_H3K56ac  <- cor.test(All_zdata_new[,14], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H3K56ac_H4K12ac   <- cor.test(All_zdata_new[,14], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x14 <- c(H3K56ac_H3K4me1$estimate, H3K56ac_H3K4me2$estimate, H3K56ac_H3K4me3$estimate, H3K56ac_H3K36me1$estimate,H3K56ac_H3K36me2$estimate, H3K56ac_H3K36me3$estimate, 
        H3K56ac_H3K9me1$estimate, H3K56ac_H3K9me2$estimate, H3K56ac_H3K9me3$estimate,H3K56ac_H3K27me3$estimate,
        H3K56ac_H3K4ac$estimate, H3K56ac_H3K9ac$estimate, H3K56ac_H3K27ac$estimate, H3K56ac_H3K56ac$estimate, H3K56ac_H4K12ac$estimate)


H4K12ac_H3K4me1  <- cor.test(All_zdata_new[,15], All_zdata_new[,1],alternative = "two.sided", method = "spearman"); 
H4K12ac_H3K4me2  <- cor.test(All_zdata_new[,15], All_zdata_new[,2],alternative = "two.sided", method = "spearman")
H4K12ac_H3K4me3  <- cor.test(All_zdata_new[,15], All_zdata_new[,3],alternative = "two.sided", method = "spearman")
H4K12ac_H3K36me1 <- cor.test(All_zdata_new[,15], All_zdata_new[,4],alternative = "two.sided", method = "spearman")
H4K12ac_H3K36me2 <- cor.test(All_zdata_new[,15], All_zdata_new[,5],alternative = "two.sided", method = "spearman")
H4K12ac_H3K36me3 <- cor.test(All_zdata_new[,15], All_zdata_new[,6],alternative = "two.sided", method = "spearman")
H4K12ac_H3K9me1   <- cor.test(All_zdata_new[,15], All_zdata_new[,7],alternative = "two.sided", method = "spearman")
H4K12ac_H3K9me2   <- cor.test(All_zdata_new[,15], All_zdata_new[,8],alternative = "two.sided", method = "spearman")
H4K12ac_H3K9me3   <- cor.test(All_zdata_new[,15], All_zdata_new[,9],alternative = "two.sided", method = "spearman")
H4K12ac_H3K27me3   <- cor.test(All_zdata_new[,15], All_zdata_new[,10],alternative = "two.sided", method = "spearman")
H4K12ac_H3K4ac   <- cor.test(All_zdata_new[,15], All_zdata_new[,11],alternative = "two.sided", method = "spearman")
H4K12ac_H3K9ac   <- cor.test(All_zdata_new[,15], All_zdata_new[,12],alternative = "two.sided", method = "spearman")
H4K12ac_H3K27ac  <- cor.test(All_zdata_new[,15], All_zdata_new[,13],alternative = "two.sided", method = "spearman")
H4K12ac_H3K56ac  <- cor.test(All_zdata_new[,15], All_zdata_new[,14],alternative = "two.sided", method = "spearman")
H4K12ac_H4K12ac   <- cor.test(All_zdata_new[,15], All_zdata_new[,15],alternative = "two.sided", method = "spearman")
x15 <- c(H4K12ac_H3K4me1$estimate, H4K12ac_H3K4me2$estimate, H4K12ac_H3K4me3$estimate, H4K12ac_H3K36me1$estimate,H4K12ac_H3K36me2$estimate, H4K12ac_H3K36me3$estimate, 
        H4K12ac_H3K9me1$estimate, H4K12ac_H3K9me2$estimate, H4K12ac_H3K9me3$estimate,H4K12ac_H3K27me3$estimate,
        H4K12ac_H3K4ac$estimate, H4K12ac_H3K9ac$estimate, H4K12ac_H3K27ac$estimate, H4K12ac_H3K56ac$estimate, H4K12ac_H4K12ac$estimate)

data <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)
colnames(data) <- c("H3K4me1","H3K4me2","H3K4me3","H3K36me1","H3K36me2","H3K36me3","H3K9me1","H3K9me2","H3K9me3","H3K27me3","H3K4ac","H3K9ac","H3K27ac","H3K56ac","H4K12ac")
rownames(data) <- c("H3K4me1","H3K4me2","H3K4me3","H3K36me1","H3K36me2","H3K36me3","H3K9me1","H3K9me2","H3K9me3","H3K27me3","H3K4ac","H3K9ac","H3K27ac","H3K56ac","H4K12ac")
write.table(data,"WT_H3K4me1_enriched_genes_modification_z-normalized数据计算spearman相关系数.xls",sep="\t")



all <- read.delim("WT_H3K4me1_enriched_genes_modification_z-normalized数据计算spearman相关系数.xls",header = T)
library(ComplexHeatmap)  #提供Heatmap函数
library(circlize)  #提供colorRamp2函数
library(RColorBrewer) #提供brewer.pal函数
heatmatrix <- as.matrix(all)
#使用cell_fun在热图显示数值
#使用cluster_rows = F, cluster_columns = F不进行聚类，热图顺序就会按照输入文件顺序,使用,width = unit(3, "cm")改变文件宽度
Heatmap(heatmatrix, name = "cor",col=colorRamp2(c(-1, 0, 1), c("blue", "#EEEEEE", "red")),cluster_rows = F, cluster_columns = F, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", heatmatrix[i, j]), x, y, gp = gpar(fontsize = 8))})
#输出大小为8*8inches


