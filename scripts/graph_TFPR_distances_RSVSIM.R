#!/usr/bin/R
## this script plots the TPR and FDR for each the telocate  and temp results of RSVSIM
## USE: (navigate to directory with BEDCOMAPRE results files) graph_TFPN_distances_RSVSIM.R
#install.packages("tidyr")
#install.packages("cowplot")
#library(cowplot)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)
directory = getwd()
setwd(directory)
##REMOVE BELOW LATER
setwd("/Users/kristen/Desktop/Fig_p")
summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)
#print(summarydata)

## CHECK NAME ORDER!!!!!!!
colnames(summarydata)[3] <- "Read_Support"

#######################################################################################

#                                 TEMP

#######################################################################################

FAM <- summarydata[ summarydata$fam=="family_aware", ]
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD<- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F_error")& summarydata$fam=="family_aware",]


merged<-merge(TEMP,SD, by="Read_Support")
#collapse 2 columns to prep for faceting
merged<-merged %>%
  gather(rate_name, rate_value, TPR.x:FDR.x)
#separate TPR and FDR and associate with the proper error number
TPR<- merged[(merged$rate_name=="TPR.x"),c(1:7,8,10:11)]
FDR<- merged[(merged$rate_name=="FDR.x"),c(1:7,9:11)]
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"
final <-merge(TPR, FDR, all=TRUE)

method_names <- list(
  'TPR.x'="TPR",
  'FDR.x'="FDR"
)

method_labeller <- function(variable,value){
  if (variable=='rate_name') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}
###
####

ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
library(ggsave)
##
#
a <- ggplot(data = final, aes(x = Read_Support, y=rate_value))
a <- a + geom_line(color="darkorange")+ xlim(0,50)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y",labeller=method_labeller) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
        
  #legend=FALSE +
  labs(x = "Read Support", y="") 
#a
#flip y-axis facet labels
g <- ggplotGrob(a)
g$layout[g$layout$name == "strip-right",c("l", "r")] <- 2
#grid.newpage()
grid.draw(g)
g
ggsave(g,filename="FAMILY_AWARE_TEMP_TPR_and_FDR.tiff",g,dpi=300, width=7.5,height=3.5,units="in")




#######################################################################################

#                                 TELOCATE

#######################################################################################

TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD<- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F_error")& summarydata$fam=="family_aware",]


merged<-merge(TELOCATE,SD, by="Read_Support")
#collapse 2 columns to prep for faceting
merged<-merged %>%
  gather(rate_name, rate_value, TPR.x:FDR.x)
#separate TPR and FDR and associate with the proper error number
TPR<- merged[(merged$rate_name=="TPR.x"),c(1:7,8,10:11)]
FDR<- merged[(merged$rate_name=="FDR.x"),c(1:7,9:11)]
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"
final <-merge(TPR, FDR, all=TRUE)

method_names <- list(
  'TPR.x'="TPR",
  'FDR.x'="FDR"
)

method_labeller <- function(variable,value){
  if (variable=='rate_name') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}
a <- ggplot(data = TPR, aes(x = Read_Support, y=rate_value))
a <- a + geom_line(color="slateblue1")+ xlim(0,50)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
        #background_grid() +
        
  #legend=FALSE +

  labs(x = "Read Support", y="TPR") 

a
#ggdraw(switch_axis_position(a, axis = 'y')) 
ggsave(filename="FAMILY_AWARE_TELOCATE_TPR.tiff",dpi=300, width=7.5,height=3.5,units="in")


stop("Stopping....comment this line to run alternative code")




a <- ggplot(data = final, aes(x = Read_Support, y=rate_value))
a <- a + geom_line(color="slateblue1")+ xlim(0,50)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y",labeller=method_labeller) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  #legend=FALSE +
  labs(x = "Read Support", y="") 
a
ggsave(filename="FAMILY_AWARE_TELOCATE_TPR.tiff",dpi=300, width=7.5,height=3.5,units="in")



#TEMP TPR
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TEMP_TPR.pdf")
a <- ggplot(data = TEMP, aes(x = Read_Support, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TEMP_filtered","TEMP_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()

#TELOCATE FDR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_FDR.pdf")
a <- ggplot(data = TELOCATE, aes(x = Read_Support, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()
#TELOCATE TPR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_TPR.pdf")
a <- ggplot(data = TELOCATE, aes(x = Read_Support, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()