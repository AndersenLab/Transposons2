#!/usr/bin/R
## this script plots the TPR and FDR for the telocate results of ROUND18
## USE: (navigate to directory with BEDCOMAPRE results files) graph_TFPN_distances_ROUND18.R

library(ggplot2)
directory = getwd()
setwd(directory)
summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)
#print(summarydata)

## CHECK NAME ORDER!!!!!!!
colnames(summarydata)[3] <- "Read_Support"

#######################################################################################

#                                  FAMILY_AWARE

#######################################################################################

#TELOCATE FDR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_1_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_1_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_FDR.pdf")
a <- ggplot(data = TELOCATE, aes(x = Read_Support, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()
#TELOCATE TPR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_1_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_1_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_TPR.pdf")
a <- ggplot(data = TELOCATE, aes(x = Read_Support, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()



