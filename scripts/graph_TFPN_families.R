#!/usr/bin/R
# this script creates histograms displaying the TPR/FDR and associated standard deviations for temp, telocate, and retroseq for each individual transposon family
# USE:  graph_TFPN_families.R

library(ggplot2)
directory = getwd()
setwd(directory)
summarydata <- read.table("BEDCOMPARE_FAMILY_MEANS.txt",header=TRUE)
print(summarydata)
names(summarydata)

###ALL ARE AT A CUTOFF DISTANCE OF 20 AND ARE THE FILTERED DATA
TEMP<- summarydata[grep("temp",  summarydata$M1), ]
RETROSEQ<- summarydata[grep("retroseq",  summarydata$M1), ]
TELOCATE<- summarydata[grep("telocate",  summarydata$M1), ]

################
#     TEMP
################

pdf(file = "IND_FAMILY_TPR_TEMP.pdf")
m <- ggplot(TEMP, aes(x=TPR))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_FDR_TEMP.pdf")
m <- ggplot(TEMP, aes(x=FDR))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_TPR_error_TEMP.pdf")
m <- ggplot(TEMP, aes(x=TPR_error))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_FDR_error_TEMP.pdf")
m <- ggplot(TEMP, aes(x=FDR_error))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()
################
#     RETROSEQ
################

pdf(file = "IND_FAMILY_TPR_RETROSEQ.pdf")
m <- ggplot(RETROSEQ, aes(x=TPR))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_FDR_RETROSEQ.pdf")
m <- ggplot(RETROSEQ, aes(x=FDR))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_TPR_error_RETROSEQ.pdf")
m <- ggplot(RETROSEQ, aes(x=TPR_error))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_FDR_error_RETROSEQ.pdf")
m <- ggplot(RETROSEQ, aes(x=FDR_error))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()


################
#     TELOCATE
################

pdf(file = "IND_FAMILY_TPR_TELOCATE.pdf")
m <- ggplot(TELOCATE, aes(x=TPR))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_FDR_TELOCATE.pdf")
m <- ggplot(TELOCATE, aes(x=FDR))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_TPR_error_TELOCATE.pdf")
m <- ggplot(TELOCATE, aes(x=TPR_error))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()

pdf(file = "IND_FAMILY_FDR_error_TELOCATE.pdf")
m <- ggplot(TELOCATE, aes(x=FDR_error))
m + geom_histogram(binwidth=2)+theme(text = element_text(size=20),axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))
dev.off()
####################################################################################################
####################################################################################################
####################################################################################################

TEMP<- summarydata[grep("temp",  summarydata$M1), ]
print(TEMP)
hist(TEMP$TPR,breaks=100,xlab="TPR(%)",main="TEMP")
hist(TEMP$FDR,breaks=100,xlab="FDR(%)",main="TEMP")
hist(TEMP$TPR_error,breaks=100,xlab="TPR SD",main="TEMP")
hist(TEMP$FDR_error,breaks=100,xlab="FDR SD",main="TEMP")

RETROSEQ<- summarydata[grep("retroseq",  summarydata$M1), ]
print(RETROSEQ)
hist(RETROSEQ$TPR,breaks=100,xlab="TPR(%)",main="RETROSEQ")
hist(RETROSEQ$FDR,breaks=100,xlab="FDR(%)",main="RETROSEQ")
hist(RETROSEQ$TPR_error,breaks=100,xlab="TPR SD",main="RETROSEQ")
hist(RETROSEQ$FDR_error,breaks=100,xlab="FDR SD",main="RETROSEQ")

TELOCATE<- summarydata[grep("telocate",  summarydata$M1), ]
print(TELOCATE)
hist(TELOCATE$TPR,breaks=100,xlab="TPR(%)",main="TELOCATE")
hist(TELOCATE$FDR,breaks=100,xlab="FDR(%)",main="TELOCATE")
hist(TELOCATE$TPR_error,breaks=100,xlab="TPR SD",main="TELOCATE")
hist(TELOCATE$FDR_error,breaks=100,xlab="FDR SD",main="TELOCATE")

