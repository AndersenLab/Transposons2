library(dplyr)
library(ggplot2)
library(data.table)

sessionInfo()

setwd("/Users/kristen/Desktop/most_recent")
#load('SignificantMappings.Rda')
load('SignificantMappings_Results_Activity.Rda')

transposon <- stringr::str_split_fixed(Mappings$traits, "_TRANS_",2)[,2]
Mappings$family <- transposon
caller <- stringr::str_split_fixed(Mappings$traits, "_TRANS_",2)[,1]
Mappings$method <- caller

#  BY METHOD
Mappings %>%
  filter(-log10(ps) > -log10(.05/8000))%>%
  ggplot(.)+
  aes(x=pos/1e6, fill=method)+
  geom_histogram()+
  facet_grid(.~chr,scale="free_x")+
  labs(x="Chromosome Position (Mb)", y="QTL Count")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        legend.key.size = unit(.2, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=9))+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="QTL_all_method.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")


# 3X BY METHOD
Mappings %>%
  filter(-log10(ps) > -log10(.05/8000))%>%
  ggplot(.)+
  aes(x=pos/1e6, fill=method)+
  geom_histogram()+
  facet_grid(method~chr,scale="free_x")+
  labs(x="Chromosome Position (Mb)", y="QTL Count")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        legend.text=element_text(size=9))+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="QTL_all_method_3x.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

#BY FAMILY
Mappings %>%
  filter(-log10(ps) > -log10(.05/8000))%>%
  ggplot(.)+
  aes(x=pos/1e6, fill=family)+
  geom_histogram()+
  facet_grid(.~chr,scale="free_x")+
  labs(x="Chromosome Position (Mb)", y="QTL Count")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        legend.text=element_text(size=9))+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="QTL_all_family.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")


#BY BASE TRAITS ONLY (not activity)
base_traits <-Mappings[(Mappings$method=="absent"| Mappings$method=="new" |Mappings$method=="reference"), ]
base_traits %>%
  filter(-log10(ps) > -log10(.05/8000))%>%
  ggplot(.)+
  aes(x=pos/1e6, fill=method)+
  geom_histogram()+
  facet_grid(.~chr,scale="free_x")+
  labs(x="Chromosome Position (Mb)", y="QTL Count")+
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
ggsave(filename="QTL_base_traits.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

## 3X BY BASE TRAITS ONLY (not activity)

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

base_traits <-Mappings[(Mappings$method=="absent"| Mappings$method=="new" |Mappings$method=="reference"), ]
base_traits %>%
  filter(-log10(ps) > -log10(.05/8000))%>%
  ggplot(.)+
  aes(x=pos/1e6, fill=method)+
  geom_histogram()+
  facet_grid(method~chr,scale="free",labeller=method_labeller)+
  labs(x="Chromosome Position (Mb)", y="QTL Count")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        legend.position=('none'))+
  scale_fill_manual(name="",
                    labels=c("absence","insertion","reference"), 
                    values = c("darkorange", "turquoise3", "slateblue1")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
ggsave(filename="QTL_base_traits_3x.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")


###EDIT BELOW
###FILTER FOR TRAIT OF INTEREST
#filter(flights, month == 1, day == 1)





##########
#GWAS MAPPINGS
##########
selection<-filter(Mappings, -log10(ps) > -log10(.05/8000))
#for trait in selection:
for (i in unique(selection$traits.i.)){
  print(i)
  Mappings %>%
    filter(traits.i. == i)%>%
    ggplot(.)+
    aes(x=pos/1e6,y=-log10(ps),fill=method)+
    geom_point(aes( color=ifelse(-log10(ps)> -log10(.05/8000), 'red', 'black')),size=1.25)+
    facet_grid(.~chr,scale="free_x")+scale_color_identity()+
    geom_hline(aes(yintercept=-log10(.05/8000)),color="red",linetype="dashed")+
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 9, colour = "black",face="bold"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.margin = unit(.6, "lines"),
          axis.ticks =element_line(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.title=element_text(size=9),
          legend.position=('none'))+
    labs(x = "Chromosome Position (Mb)",y="-log10(p)") +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
  out_tiff <-paste(i, "tiff", sep =".")
  ggsave(filename=out_tiff,
         dpi=300,
         width=7.5,
         height=3.5,
         units="in")
}

###REMOVE LATER
Mappings %>%
  filter(traits.i. == "new_TRANS_CELE14A")%>%
  ggplot(.)+
  aes(x=pos/1e6,y=-log10(ps),fill=method)+
  geom_point(aes( color=ifelse(-log10(ps)> -log10(.05/8000), 'red', 'black')),size=1.25)+
  facet_grid(.~chr,scale="free_x")+scale_color_identity()+
  geom_hline(aes(yintercept=-log10(.05/8000)),color="red",linetype="dashed")+
  theme(strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(size = 9, colour = "black",face="bold"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
      panel.margin = unit(.6, "lines"),
      axis.ticks =element_line(colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.text.x = element_text(colour = "black"),
      axis.title=element_text(size=9),
      legend.position=('none'))+
  labs(x = "Chromosome Position (Mb)",y="-log10(p)") +
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand = c(0,0))
  out_tiff <-paste("new_TRANS_CELE14A", "tiff", sep =".")
  print(plot$ymax)
ggsave(filename=out_tiff,
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")



sigs<-filter(Mappings, -log10(ps) > -log10(.05/8000))
sigs2<-filter(sigs, chr=="IV")
print (sigs2$pos)/1e6
print(sigs2$trait)
print(sigs2)
sigs1<-filter(sigs, chr=="I")
sigs2<-filter(sigs, chr=="II")
sigs3<-filter(sigs, chr=="III")
sigs4<-filter(sigs, chr=="IV")
sigs5<-filter(sigs, chr=="V")
sigs6<-filter(sigs, chr=="X")

nonredundant1 <- unique(sigs1$pos)
nonredundant2 <- unique(sigs2$pos)
nonredundant3 <- unique(sigs3$pos)
nonredundant4 <- unique(sigs4$pos)
nonredundant5 <- unique(sigs5$pos)
nonredundant6 <- unique(sigs6$pos)
print(nonredundant1)
print(nonredundant2)
print(nonredundant3)
print(nonredundant4)
print(nonredundant5)
print(nonredundant6)

traits4 <- unique(sigs4$trait)
print(traits4)
# 33 traits and 17 SNPs
print(sigs4$trait)
print(sigs4)


