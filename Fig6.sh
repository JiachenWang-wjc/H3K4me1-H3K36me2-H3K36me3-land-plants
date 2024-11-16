#Figure 6B
#Violin-Box plot
#Box data from bw matrix.
library(tidyverse)
datag <- read.delim('box.txt',header = T)
datag <- gather(datag) %>% rename('Sample'='key','Ehd3SDG724'='value')
clean_datag <- na.omit(datag)
ggplot(clean_datag,aes(`Sample`,`Ehd3SDG724`))+ geom_violin(width=0.8,cex=0.01,trim = TRUE,color="black",aes(fill=Sample))+geom_boxplot(width=0.2,cex=0.3,outlier.shape=NA,position=position_dodge(0.5),aes(fill=Sample))+ylim(-50,50)+
scale_fill_manual(values = c("#54FF9F","#FFD700","#FF6A6A","#54FF9F","#FFD700","#FF6A6A","#54FF9F","#FFD700","#FF6A6A","#54FF9F","#FFD700","#FF6A6A"))+
#scale_fill_manual(values = rainbow(18))+
theme_bw()+theme(panel.background=element_rect(colour="black",size=1,fill="white"),axis.line = element_line(colour = "black",size=0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), 
        axis.text.x=element_text(colour="black",face="bold",size=14,angle = 45,vjust = 0.5,hjust = 0.5), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=14), 
        axis.title.x=element_text(size=15,face="bold"),
        axis.title.y=element_text(size=15,face="bold")) +
  labs(x = "Group", y = "Normalized Level of Location")     
#Save as PDF file.