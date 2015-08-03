#!/usr/bin/R
# this script plots the distrubition of transposons (totals across all samples) according to their chromosome positions
# USE: TE_density.R

library(ggplot2)
library(grid)
directory = getwd()
#setwd(directory)
###REMOVE LATER
setwd("/Users/kristen/Desktop/Fig_p")
summarydata <- read.table("all_nonredundant.txt",header=TRUE)
#summarydata <- read.table("CB4856_temp_insertion_nonredundant.bed",header=TRUE)
names(summarydata)
names(summarydata)<-c("chr","start","end","TE","support","orientation","method")

m <- ggplot(summarydata, aes(x=start/1e6,fill=method))
m + geom_bar(binwidth=1)+ 
  facet_grid(. ~ chr,scale="free_x")+
  labs(x="Chromosome Position (Mb)", y="Number of Transposons")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        legend.text=element_text(size=9))+
  scale_fill_manual(name="",
                    labels=c("absence","insertion","reference"), 
                    values = c("darkorange", "turquoise3", "slateblue1")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="chromosome_distribution.tiff",
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
  if (variable=='method') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}

m <- ggplot(summarydata, aes(x=start/1e6,fill=method))
m + geom_bar(binwidth=1)+ 
  facet_grid(method ~ chr,scale="free_x",labeller=method_labeller)+
  labs(x="Chromosome Position (Mb)", y="Number of Transposons")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        legend.text=element_text(size=9))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1")) +
  guides(fill=FALSE) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="chromosome_distribution_3x.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

##same y axis scale?