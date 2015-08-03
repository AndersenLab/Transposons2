#!/usr/bin/R
# this script
# 1) plots a stacked histogram representing the total number of transposons found in each strain per insertions, references, and absences
# 2) plots individual histogram representing the total number of transposons found in each strain per insertions, references, and absences
# 3) plots total transposons vs strain per insertions, references, and absences
# USE: Te_totals_distribution.R


#HISTOGRAMS
library(ggplot2)
directory = getwd()
setwd(directory)
##REMOVE BELOW LATER
setwd("/Users/kristen/Desktop/Fig_p")
summarydata <- read.table("FINAL_RESULTS.txt",header=TRUE)



#STACKED
m <- ggplot(summarydata, aes(x=total_tes,fill=method))
m + geom_bar(binwidth=20) +
theme(legend.position=c(.90,0.75),
      legend.background = element_rect(fill=FALSE),
      legend.text=element_text(size=9),
      panel.background = element_rect(fill = "white"),
      axis.ticks =element_line(colour = "black"),
      axis.text.y = element_text(colour = "black",size=9),
      axis.text.x = element_text(colour = "black",size=9),
      axis.line=element_line(linetype="solid"),
      axis.title=element_text(size=9))+
scale_fill_manual(name="",
                  labels=c("absence","insertion","reference"), 
                  values = c("darkorange", "turquoise3", "slateblue1")) +
labs(x = "Transposons Per Strain", y = "Count") +
scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="Histogram_of_Number_of_Transposons_in_a_Strain.tiff",
       dpi=300, 
       width=7.5,
       height=3.5,
       units="in")

#3X
method_names <- list(
  'absent'="absence",
  'new'="insertion",
  'reference'="reference"
)

method_labeller <- function(variable,value){
  return(method_names[value])
}

m <- ggplot(summarydata, aes(x=total_tes),fill=method))
m + geom_bar(binwidth=20) +
  facet_grid(. ~ method,labeller=method_labeller )+
  theme(strip.text.x = element_text(size = 12, colour = "black"),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1")) +
  guides(fill=FALSE) +
  labs(x = "Transposons Per Strain", y = "Count") +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="Histogram_of_Number_of_Transposons_in_a_Strain_3x.tiff",
       dpi=300, 
       width=7.5,
       height=3.5,
       units="in")

#########################
#INDIVIDUAL INSERTION,ABSENCE,REFERENCE
pdf(file = "insertions.pdf")
insertions<-summarydata[summarydata$method=="new",]
m <- ggplot(insertions, aes(x=total_tes))
m + geom_bar(binwidth=25)
dev.off()

absences<-summarydata[summarydata$method=="absent",]
pdf(file = "absences.pdf")
m <- ggplot(absences, aes(x=total_tes))
m + geom_bar(binwidth=25)
dev.off()

references<-summarydata[summarydata$method=="reference",]
pdf(file = "references.pdf")
m <- ggplot(references, aes(x=total_tes))
m + geom_bar(binwidth=25)
dev.off()

#TRANSPOSONS vs STRAINS
names(summarydata)
#INSERTIONS
insertions<-summarydata[summarydata$method=="new",]
insertions<-(insertions[ order(insertions$total_tes), ])
plot(insertions$total_tes~insertions$sample)
pdf(file = "insertions2.pdf")
m <- ggplot(insertions, aes(x=reorder(insertions$sample,insertions$total_tes), y=insertions$total_tes)) 
m + geom_line() + stat_smooth(colour='blue')+aes(group=1)
dev.off()
#ABSENCES
absences<-summarydata[summarydata$method=="absent",]
absences<-(absences[ order(absences$total_tes), ])
plot(absences$total_tes~absences$sample)
pdf(file = "absences2.pdf")
m <- ggplot(absences, aes(x=reorder(absences$sample,absences$total_tes), y=absences$total_tes)) 
m + geom_line() + stat_smooth(colour='blue')+aes(group=1)
dev.off()
#REFERENCE
references<-summarydata[summarydata$method=="reference",]
references<-(references[ order(references$total_tes), ])
plot(references$total_tes~references$sample)
pdf(file = "references2.pdf")
m <- ggplot(references, aes(x=reorder(references$sample,references$total_tes), y=references$total_tes)) 
m + geom_line() + stat_smooth(colour='blue')+aes(group=1)
dev.off()