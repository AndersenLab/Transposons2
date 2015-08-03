#!/usr/bin/R
## this script plots the TPR and FDR for each transposon detection method
## produces graphs for both overall transposon detection and family-aware detection
## plots each method seaparately and on the same graph
## creates histograms displaying the distance between the position of the simulated TE and the position of the detected transposon
## NOTE: set x limits separately based on histograms

## USE: (navigate to directory with BEDCOMAPRE results files) graph_TFPN_distances.R

library(ggplot2)
directory = getwd()
setwd(directory)
summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)
#print(summarydata)

## CHECK NAME ORDER!!!!!!!
#######################################################################################

#                                  OVERALL

#######################################################################################
OVERALL <- summarydata[ summarydata$fam=="overall", ]
#TEMP FDR
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_NF" )& summarydata$fam=="overall", ]
pdf(file = "OVERALL_TEMP_FDR.pdf")
a <- ggplot(data = TEMP, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TEMP_filtered","TEMP_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#TEMP TPR
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_NF" )& summarydata$fam=="overall", ]
pdf(file = "OVERALL_TEMP_TPR.pdf")
a <- ggplot(data = TEMP, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TEMP_filtered","TEMP_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#TELOCATE FDR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="overall", ]
pdf(file = "OVERALL_TELOCATE_FDR.pdf")
a <- ggplot(data = TELOCATE, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#TELOCATE TPR
TELCOATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telcoate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="overall", ]
pdf(file = "OVERALL_TELOCATE_TPR.pdf")
a <- ggplot(data = TELOCATE, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#RETROSEQ FDR
RETROSEQ <- summarydata[(summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_NF" )& summarydata$fam=="overall", ]
pdf(file = "OVERALL_RETROSEQ_FDR.pdf")
a <- ggplot(data = RETROSEQ, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#RETROSEQ TPR
RETROSEQ <- summarydata[(summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_NF" )& summarydata$fam=="overall", ]
pdf(file = "OVERALL_RETROSEQ_TPR.pdf")
a <- ggplot(data = RETROSEQ, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()
#######################################################################################

#                                  FAMILY_AWARE

#######################################################################################
FAM <- summarydata[ summarydata$fam=="family_aware", ]
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TEMP_FDR.pdf")
a <- ggplot(data = TEMP, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TEMP_filtered","TEMP_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#TEMP TPR
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TEMP_TPR.pdf")
a <- ggplot(data = TEMP, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,100)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TEMP_filtered","TEMP_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#TELOCATE FDR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_FDR.pdf")
a <- ggplot(data = TELOCATE, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,400)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()
#TELOCATE TPR
TELCOATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telcoate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_TPR.pdf")
a <- ggplot(data = TELOCATE, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,400)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#RETROSEQ FDR
RETROSEQ <- summarydata[(summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_RETROSEQ_FDR.pdf")
a <- ggplot(data = RETROSEQ, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()

#RETROSEQ TPR
RETROSEQ <- summarydata[(summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_RETROSEQ_TPR.pdf")
a <- ggplot(data = RETROSEQ, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_unfiltered"), values = c("red", "darkturquoise"))
a
dev.off()
#######################################################################################

#                                   ALL

#######################################################################################
#OVERALL FDR
N2_OVERALL <- OVERALL[(grep("N2", OVERALL$M2)), ]
pdf(file = "OVERALL_ALL_FDR.pdf")
a <- ggplot(data = N2_OVERALL, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,400)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_filtered_error","RETROSEQ_unfiltered","RETROSEQ_unfiltered_error", "TELOCATE_filtered","TELOCATE_filtered_error","TELOCATE_unfiltered","TELCOATE_unfilted_error","TEMP_filtered","TEMP_filtered_error", "TEMP_unfiltered","TEMP_unfiltered_error"), values = c("red", "orange","green","lightgreen","blue","purple","red2", "orange2","green2","red3","blue2","purple2"))
a
dev.off()

#OVERALL TPR
N2_OVERALL <- OVERALL[(grep("N2", OVERALL$M2)), ]
pdf(file = "OVERALL_ALL_TPR.pdf")
a <- ggplot(data = N2_OVERALL, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,400)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_filtered_error","RETROSEQ_unfiltered","RETROSEQ_unfiltered_error", "TELOCATE_filtered","TELOCATE_filtered_error","TELOCATE_unfiltered","TELCOATE_unfilted_error","TEMP_filtered","TEMP_filtered_error", "TEMP_unfiltered","TEMP_unfiltered_error"), values = c("red", "orange","green","lightgreen","blue","purple","red2", "orange2","green2","red3","blue2","purple2"))
a
dev.off()

#FAMILY_AWARE FDR
FAM <- summarydata[ summarydata$fam=="family_aware", ]
N2_FAM <- FAM[(grep("N2", FAM$M2)), ]
pdf(file = "FAMILY_AWARE_ALL_FDR.pdf")
a <- ggplot(data = N2_FAM, aes(x = DistanceCutoff, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,400)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_filtered_error","RETROSEQ_unfiltered","RETROSEQ_unfiltered_error", "TELOCATE_filtered","TELOCATE_filtered_error","TELOCATE_unfiltered","TELCOATE_unfilted_error","TEMP_filtered","TEMP_filtered_error", "TEMP_unfiltered","TEMP_unfiltered_error"), values = c("red", "orange","green","lightgreen","blue","purple","red2", "orange2","green2","red3","blue2","purple2"))
a
dev.off()

#FAMILY_AWARE_TPR
N2_FAM <- FAM[(grep("N2", FAM$M2)), ]
pdf(file = "FAMILY_AWARE_ALL_TPR.pdf")
a <- ggplot(data = N2_FAM, aes(x = DistanceCutoff, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,400)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("RETROSEQ_filtered","RETROSEQ_filtered_error","RETROSEQ_unfiltered","RETROSEQ_unfiltered_error", "TELOCATE_filtered","TELOCATE_filtered_error","TELOCATE_unfiltered","TELCOATE_unfilted_error","TEMP_filtered","TEMP_filtered_error", "TEMP_unfiltered","TEMP_unfiltered_error"), values = c("red", "orange","green","lightgreen","blue","purple","red2", "orange2","green2","red3","blue2","purple2"))
a
dev.off()


#######################################################################################

#                                   HISTOGRAMS

#######################################################################################

#         TEST RUN BELOW!!!

directory = getwd()
setwd(directory)
distance_data<-read.table("SCs_ALL_filter")
print(distance_data)
names(distance_data)
colnames(distance_data) <- c("chr1","start1","nd1","family","null1","strand1","chr2","start2","end2","method","null2","strand2","distance")

TEMP<- distance_data[grep("temp",  distance_data$method), ]
pdf(file = "Distance_Histogram_TEMP_F.pdf")
m <- ggplot(TEMP, aes(x=distance))
#m + geom_histogram()
m + geom_histogram(binwidth=2)+scale_x_continuous(limits = c(0, 100))+ggtitle("TEMP")+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

RETROSEQ<- distance_data[grep("retroseq",  distance_data$method), ]
pdf(file = "Distance_Histogram_RETROSEQ_F.pdf")
m <- ggplot(RETROSEQ, aes(x=distance))
#m + geom_histogram()
m + geom_histogram(binwidth=2)+scale_x_continuous(limits = c(0, 50))+ggtitle("RETROSEQ")+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

TELOCATE<- distance_data[grep("telocate",  distance_data$method), ]
pdf(file = "Distance_Histogram_TELCOATE_F.pdf")
m <- ggplot(TELOCATE, aes(x=distance))
#m + geom_histogram()
m + geom_histogram(binwidth=2)+scale_x_continuous(limits = c(0, 400))+ggtitle("TELOCATE")+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()




distance_data<-read.table("SCs_ALL_non_filter")
print(distance_data)
names(distance_data)
colnames(distance_data) <- c("chr1","start1","nd1","family","null1","strand1","chr2","start2","end2","method","null2","strand2","distance")

TEMP<- distance_data[grep("temp",  distance_data$method), ]
pdf(file = "Distance_Histogram_TEMP_NF.pdf")
m <- ggplot(TEMP, aes(x=distance))
m + geom_histogram()
m + geom_histogram(binwidth=2)+scale_x_continuous(limits = c(0, 100))+ggtitle("TEMP")+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

RETROSEQ<- distance_data[grep("retroseq",  distance_data$method), ]
pdf(file = "Distance_Histogram_RETROSEQ_NF.pdf")
m <- ggplot(RETROSEQ, aes(x=distance))
m + geom_histogram()
m + geom_histogram(binwidth=2)+scale_x_continuous(limits = c(0, 50))+ggtitle("RETROSEQ")+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

TELOCATE<- distance_data[grep("telocate",  distance_data$method), ]
pdf(file = "Distance_Histogram_TELCOATE_NF.pdf")
m <- ggplot(TELOCATE, aes(x=distance))
m + geom_histogram()
m + geom_histogram(binwidth=2)+scale_x_continuous(limits = c(0, 400))+ggtitle("TELOCATE")+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()