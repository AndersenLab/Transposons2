#!/usr/bin/R
## this script plots total TE counts vs coverage
## USE: coverage.R

setwd("/Users/kristen/Desktop/Fig_p")
summarydata <- read.table("coverage_and_te_counts.txt",header=TRUE)
library(tidyr)
library(ggplot2)
library(dplyr)

print(summarydata)
with(summarydata, cor.test(absence, coverage))
with(summarydata, cor.test(insertion, coverage))
with(summarydata, cor.test(reference, coverage))
with(summarydata, cor(absence, coverage))
with(summarydata, cor(insertion, coverage))
with(summarydata, cor(reference, coverage))

fit<-lm(absence ~ coverage, data=summarydata)
plot(summarydata$coverage,summarydata$absence,ylim=c(0,300))
abline(fit)
summary(fit)
fit<-lm(insertion ~ coverage, data=summarydata)
plot(summarydata$coverage,summarydata$insertion)
abline(fit)
summary(fit)
fit<-lm(reference ~ coverage, data=summarydata)
plot(summarydata$coverage,summarydata$reference)
abline(fit)
summary(fit)

summary(test)

plot(summarydata$absence~summarydata$coverage)
plot(summarydata$insertion~summarydata$coverage)
plot(summarydata$reference~summarydata$coverage)




fit<-lm(absence ~ coverage, data=summarydata)
summary(fit)
fit<-lm(insertion ~ coverage, data=summarydata)
summary(fit)
fit<-lm(reference ~ coverage, data=summarydata)
summary(fit)

names(data)
data<-summarydata %>%
  gather(method, total_tes, absence:insertion:reference)
names(data)
print(data)

#LINE
a <- ggplot(data = data, aes(x = coverage,y=total_tes, colour=method))
a <- a + geom_line(size=1)+
  theme(legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_color_manual(name="",
                    #labels=c("absence","insertion","reference"),
                    values = c("darkorange", "turquoise3", "slateblue1"))+
                    labs(x = "Depth of Coverage", y="Number of Transposons")
a
ggsave(filename="totals_vs_coverage_line.tiff",dpi=300, width=7.5,height=3.5,units="in")




#3X LINE
a <- ggplot(data = data, aes(x = coverage,y=total_tes, colour=method))
a <- a + geom_line(size=1)+
  facet_grid(method ~ .,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        legend.position=('none'),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_color_manual(name="",
                     #labels=c("absence","insertion","reference"),
                     values = c("darkorange", "turquoise3", "slateblue1"))+
  labs(x = "Depth of Coverage", y="Number of Transposons")
a
ggsave(filename="totals_vs_coverage_line_3x.tiff",dpi=300, width=7.5,height=3.5,units="in")


#SCATTER
a <- ggplot(data = data, aes(x = coverage,y=total_tes, colour=method))
a <- a + geom_point(size=1)+
  geom_smooth(method=lm,se=FALSE)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=c(.90,0.85),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9))+
  scale_color_manual(name="",
                     #labels=c("absence","insertion","reference"),
                     values = c("darkorange", "turquoise3", "slateblue1"))+
  labs(x = "Depth of Coverage", y="Number of Transposons")
a
ggsave(filename="totals_vs_coverage.tiff",dpi=300, width=7.5,height=3.5,units="in")


#3X
a <- ggplot(data = data, aes(x = coverage,y=total_tes, colour=method))
a <- a + geom_point(size=1)+
  geom_smooth(method=lm,se=FALSE)+
  facet_grid(method ~ .,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        legend.position=('none'),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_color_manual(name="",
                     values = c("darkorange", "turquoise3", "slateblue1"))+
  labs(x = "Depth of Coverage", y="Number of Transposons")
a
ggsave(filename="totals_vs_coverage_3x.tiff",dpi=300, width=7.5,height=3.5,units="in")



###PDF
a <- ggplot(data = data, aes(x = coverage,y=total_tes, colour=method))
a <- a + geom_point(size=1)+
  geom_smooth(method=lm,se=FALSE)+
  facet_grid(method ~ .,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        legend.position=('none'),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_color_manual(name="",
                     #labels=c("absence","insertion","reference"),
                     values = c("darkorange", "turquoise3", "slateblue1"))+
  labs(x = "Depth of Coverage", y="Number of Transposons")
a
ggsave(filename="totals_vs_coverage_3x.pdf",dpi=300, width=7.5,height=3.5,units="in")