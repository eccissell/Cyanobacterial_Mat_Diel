library(dplyr)
library(ggplot2)
library(viridis)
library(tidyverse)
library(gdata)
library(reshape2)
library(vegan)
library(conflicted)
library(ggpubr)
library(lme4)
library(lmerTest)
library(gridExtra)
library(nlme)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(ggsignif)
library(car)
conflict_prefer("summarise",'dplyr')
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")

#####DIEL RNA and DNA DATA ALIGNED TO SECONDARY METABOLITE REGIONS IN DOMINANT CYANO MAG
#Read in all the counts (TPM)
SEC_D1=read_tsv('/Users/cissell/Desktop/D1_sec_metab_tpm.tsv')
SEC_D1=as.data.frame(SEC_D1)
SEC_D2=read_tsv('/Users/cissell/Desktop/D2_sec_metab_tpm.tsv')
SEC_D2=as.data.frame(SEC_D2)
SEC_D3=read_tsv('/Users/cissell/Desktop/D3_sec_metab_tpm.tsv')
SEC_D3=as.data.frame(SEC_D3)
SEC_D4=read_tsv('/Users/cissell/Desktop/D4_sec_metab_tpm.tsv')
SEC_D4=as.data.frame(SEC_D4)

SEC_R1=read_tsv('/Users/cissell/Desktop/R1_sec_metab_tpm.tsv')
SEC_R1=as.data.frame(SEC_R1)
SEC_R2=read_tsv('/Users/cissell/Desktop/R2_sec_metab_tpm.tsv')
SEC_R2=as.data.frame(SEC_R2)
SEC_R3=read_tsv('/Users/cissell/Desktop/R3_sec_metab_tpm.tsv')
SEC_R3=as.data.frame(SEC_R3)
SEC_R4=read_tsv('/Users/cissell/Desktop/R4_sec_metab_tpm.tsv')
SEC_R4=as.data.frame(SEC_R4)

#Build master DNA
LSEC_D1=gather(SEC_D1,Source,TPM,D1_1:D1_5,factor_key=T)
LSEC_D2=gather(SEC_D2,Source,TPM,D2_1:D2_5,factor_key=T)
LSEC_D3=gather(SEC_D3,Source,TPM,D3_1:D3_5,factor_key=T)
LSEC_D4=gather(SEC_D4,Source,TPM,D4_1:D4_5,factor_key=T)

#Build master RNA
LSEC_R1=gather(SEC_R1,Source,TPM,R1_1:R1_5,factor_key=T)
LSEC_R2=gather(SEC_R2,Source,TPM,R2_1:R2_5,factor_key=T)
LSEC_R3=gather(SEC_R3,Source,TPM,R3_1:R3_5,factor_key=T)
LSEC_R4=gather(SEC_R4,Source,TPM,R4_1:R4_5,factor_key=T)

#Build Master
SEC_D=rbind(LSEC_D1,LSEC_D2,LSEC_D3,LSEC_D4)
SEC_R=rbind(LSEC_R1,LSEC_R2,LSEC_R3,LSEC_R4)
SEC=rbind(SEC_D,SEC_R)

#Add new columnns for sorting
SEC=SEC %>%
  mutate(Molecule=case_when(
    str_detect(Source,"R") ~ "R",
    str_detect(Source,"D") ~ "D",
  ))
SEC=SEC %>%
  mutate(Mat=case_when(
    str_detect(Source,"R1") ~ "1",
    str_detect(Source,"D1") ~ "1",
    str_detect(Source,"D2") ~ "2",
    str_detect(Source,"R2") ~ "2",
    str_detect(Source,"D3") ~ "3",
    str_detect(Source,"R3") ~ "3",
    str_detect(Source,"D4") ~ "4",
    str_detect(Source,"R4") ~ "4",
  ))
SEC$Source=as.character(SEC$Source)
SEC=SEC %>%
  mutate(Time=case_when(
    endsWith(Source,"1") ~ "1",
    endsWith(Source,"2") ~ "2",
    endsWith(Source,"3") ~ "3",
    endsWith(Source,"4") ~ "4",
    endsWith(Source,"5") ~ "5",
  ))
SEC$Time=as.factor(SEC$Time)
SEC$Mat=as.factor(SEC$Mat)

#Fit regressions
#RNA subset to NRPS
sec_lm1=lm(TPM~factor(Time)factor(Mat),data=subset(SEC,Type=="NRPS" | Molecule=="R"))
summary(sec_lm1)
check_gaussian_model1 <- simulateResiduals(fittedModel = sec_lm1, n = 500)
plot(check_gaussian_model1)
hist(resid(sec_lm1))
#Try gamma distribution
sec_glm1=glm(TPM~factor(Time)+Mat,data=subset(SEC,Type=="NRPS" | Molecule=="R"), family=Gamma(link="log"))
summary(sec_glm1)
#See if interaction can be ignored
sec_glm1_inter=glm(TPM~factor(Time)*factor(Mat),data=subset(SEC,Type=="NRPS" | Molecule=="R"), family=Gamma(link="log"))
#Test
anova(sec_glm1,sec_glm1_inter,test="LRT")
#Yes it can be ignored
#Check assumptions
library("DHARMa")
check_gamma_model1 <- simulateResiduals(fittedModel = sec_glm1, n = 500)
plot(check_gamma_model1)
#All good
Anova(sec_glm1,type="III")
#Get pairwise contrasts
library(multcomp)
summary(glht(sec_glm1, linfct=mcp(Mat="Tukey")))


#RNA all
sec_glm2=glm(TPM~factor(Time)+Mat,data=subset(SEC,Molecule=="R"), family=Gamma(link="log"))
summary(sec_glm2)
#See if interaction can be ignored
sec_glm2_inter=glm(TPM~factor(Time)*Mat,data=subset(SEC, Molecule=="R"), family=Gamma(link="log"))
#Test
anova(sec_glm2,sec_glm2_inter,test="LRT")
#Can be ignored
#Check assumptions
library("DHARMa")
check_gamma_model2 <- simulateResiduals(fittedModel = sec_glm2, n = 500)
plot(check_gamma_model2)
#All good
Anova(sec_glm2,type="III")
#Get pairwise contrasts
summary(glht(sec_glm2, linfct=mcp(Mat="Tukey")))


#DNA subset to NRPS
sec_glm3=glm(TPM~factor(Time)+Mat,data=subset(SEC,Type=="NRPS" | Molecule=="D"), family=Gamma(link="log"))
summary(sec_glm3)
#See if interaction can be ignored
sec_glm3_inter=glm(TPM~factor(Time)*factor(Mat),data=subset(SEC,Type=="NRPS" | Molecule=="D"), family=Gamma(link="log"))
#Test
anova(sec_glm3,sec_glm3_inter,test="LRT")
#Check assumptions
library("DHARMa")
check_gamma_model3 <- simulateResiduals(fittedModel = sec_glm3, n = 500)
plot(check_gamma_model3)
#Not great but seems ok visually
Anova(sec_glm3,type="III")
#Get pairwise contrasts
summary(glht(sec_glm3, linfct=mcp(Mat="Tukey")))



#DNA all
sec_glm4=glm(TPM~factor(Time)+Mat,data=subset(SEC,Molecule=="D"), family=Gamma(link="log"))
summary(sec_glm4)
#See if interaction can be ignored
sec_glm4_inter=glm(TPM~factor(Time)*factor(Mat),data=subset(SEC, Molecule=="D"), family=Gamma(link="log"))
#Test
anova(sec_glm4,sec_glm4_inter,test="LRT")
#Check assumptions
library("DHARMa")
check_gamma_model4 <- simulateResiduals(fittedModel = sec_glm4, n = 500)
plot(check_gamma_model4)
#Not great but seems ok visually
Anova(sec_glm4,type="III")
#Get pairwise contrasts
summary(glht(sec_glm4, linfct=mcp(Mat="Tukey")))


#This is all depricated now, but supports the GLM inference so that's cool.
#Run some models to look for significant changes across time
kruskal.test(TPM~Time,data=subset(SEC,Type=="NRPS" | Molecule=="R"))
kruskal.test(TPM~Time,data=subset(SEC,Type=="NRPS" | Molecule=="D"))
kruskal.test(TPM~Time,data=subset(SEC,Molecule=="R"))
kruskal.test(TPM~Time,data=subset(SEC,Molecule=="D"))
#NO SIG DIFS
#Now try spatial dif
kruskal.test(TPM~Mat,data=subset(SEC,Molecule=="R")) #Sig
kruskal.test(TPM~Mat,data=subset(SEC,Molecule=="D")) #Sig
kruskal.test(TPM~Mat,data=subset(SEC,Type=="NRPS" | Molecule=="R")) #Sig
kruskal.test(TPM~Mat,data=subset(SEC,Type=="NRPS" | Molecule=="D")) #Sig
#Post-Hoc pairwise
#Create subset dataframes for these tests
sub_R=subset(SEC,Molecule=="R")
sub_D=subset(SEC,Molecule=="D")
sub_R_N=subset(SEC,Type=="NRPS" | Molecule=="R")
sub_D_N=subset(SEC,Type=="NRPS" | Molecule=="D")
pairwise.wilcox.test(sub_R$TPM, sub_R$Mat, p.adjust.method="bonferroni")
pairwise.wilcox.test(sub_D$TPM, sub_D$Mat, p.adjust.method="bonferroni")
pairwise.wilcox.test(sub_D_N$TPM, sub_D_N$Mat, p.adjust.method="bonferroni")
pairwise.wilcox.test(sub_R_N$TPM, sub_R_N$Mat, p.adjust.method="bonferroni")
?pairwise.wilcox.test

#Boxplots
overall_box=ggplot(data=SEC,aes(y=TPM,x=Time, fill=Molecule)) +
  geom_boxplot(outlier.shape = NA,
               size=1,
               position = position_dodge2(preserve = "single"),
               alpha=0.9) +
  geom_point(aes(fill=Molecule),position=position_jitterdodge(),colour="black",pch=21,
             alpha=0.6,
             size=1)+
  theme_classic()+
  scale_fill_manual(name = "Read Set", labels = c("DNA", "RNA"),
                    values=c("#455765","#cd9fb2"))+
  scale_color_manual(values=c("#455765","#cd9fb2"))+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_y_continuous(name='')+
  scale_x_discrete(breaks=c(1,2,3,4,5), name='',
                   labels=c("1"="09:00","2"="15:00","3"="21:00","4"="03:00","5"="09:00"))+
  theme(legend.position = "none")
overall_box

overall_box_mat=ggplot(data=SEC,aes(y=TPM,x=Mat, fill=Molecule)) +
  geom_boxplot(outlier.shape = NA,
               size=1,
               alpha=0.9) +
  geom_point(aes(fill=Molecule),position=position_jitterdodge(),colour="black",pch=21,
             alpha=0.6,
             size=1)+
  theme_classic()+
  scale_fill_manual(name = "Read Set", labels = c("DNA", "RNA"),
                    values=c("#455765","#cd9fb2"))+
  scale_color_manual(values=c("#455765","#cd9fb2"))+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_y_continuous(name='')+
  scale_x_discrete(name='')+
  geom_signif(color="#333333",
              xmin=2,
              xmax=3,
              y_position=220000,
              annotation="**",
              size=1,
              textsize=7,
              vjust=0.6) +
  geom_signif(color="#333333",
              xmin=2,
              xmax=4,
              y_position=240000,
              annotation="**",
              size=1,
              textsize=7,
              vjust=0.6)+
  theme(legend.position = "none")
overall_box_mat

NRPS_box=ggplot(data=subset(SEC,Type=="NRPS"),aes(y=TPM,x=Time, fill=Molecule)) +
  geom_boxplot(outlier.shape = NA,
               size=1,
               position = position_dodge2(preserve = "single"),
               alpha=0.9) +
  geom_point(aes(fill=Molecule),position=position_jitterdodge(),colour="black",pch=21,
             alpha=0.6,
             size=1)+
  theme_classic()+
  scale_fill_manual(name = "Read Set", labels = c("DNA", "RNA"),
                    values=c("#455765","#cd9fb2"))+
  scale_color_manual(values=c("#455765","#cd9fb2"))+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_y_continuous(name='')+
  scale_x_discrete(breaks=c(1,2,3,4,5), name="Time Point",
                   labels=c("1"="09:00","2"="15:00","3"="21:00","4"="03:00","5"="09:00"))+
  theme(legend.position = "none")
NRPS_box

NRPS_box_mat=ggplot(data=subset(SEC,Type=="NRPS"),aes(y=TPM,x=Mat, fill=Molecule)) +
  geom_boxplot(outlier.shape = NA,
               size=1,
               alpha=0.9) +
  geom_point(aes(fill=Molecule),position=position_jitterdodge(),colour="black",pch=21,
             alpha=0.6,
             size=1)+
  theme_classic()+
  scale_fill_manual(name = "Read Set", labels = c("DNA", "RNA"),
                    values=c("#455765","#cd9fb2"))+
  scale_color_manual(values=c("#455765","#cd9fb2"))+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_y_continuous(name='')+
  geom_signif(color="#cd9fb2",
              xmin=1,
              xmax=2,
              y_position=140000,
              annotation="*",
              size=1,
              textsize=7,
              vjust=0.6) +
  geom_signif(color="#cd9fb2",
              xmin=2,
              xmax=3,
              y_position=150000,
              annotation="***",
              size=1,
              textsize=7,
              vjust=0.6,
              hjust=-0.05)+
  geom_signif(color="#cd9fb2",
              xmin=2,
              xmax=4,
              y_position=160000,
              annotation="***",
              size=1,
              textsize=7,
              vjust=0.6,
              hjust=-0.5)+
  geom_signif(color="#333333",
              xmin=2,
              xmax=3,
              y_position=150000,
              annotation="**",
              size=1,
              textsize=7,
              vjust=0.6,
              hjust=1.25)+
  geom_signif(color="#333333",
              xmin=2,
              xmax=4,
              y_position=160000,
              annotation="**",
              size=1,
              textsize=7,
              vjust=0.6,
              hjust=2)+
  theme(legend.position = "none")
NRPS_box_mat

#stacked barplot of different types of secondary metabolite clusters
#Create proportions
conflict_prefer("filter", "dplyr")
sum_SEC_R=SEC%>%
  filter(Molecule=="R") %>%
  group_by(Type,Molecule) %>%
  summarise(n=sum(TPM)) %>%
  ungroup %>%
  mutate(prop=n/sum(n),
         lower = lapply(n, prop.test, n = sum(n)),
         upper = sapply(lower, function(x) x$conf.int[2]),
         lower = sapply(lower, function(x) x$conf.int[1]))

sum_SEC_D=SEC%>%
  filter(Molecule=="D") %>%
  group_by(Type,Molecule) %>%
  summarise(n=sum(TPM)) %>%
  ungroup %>%
  mutate(prop=n/sum(n),
         lower = lapply(n, prop.test, n = sum(n)),
         upper = sapply(lower, function(x) x$conf.int[2]),
         lower = sapply(lower, function(x) x$conf.int[1]))
?prop.test
#Combine into master proportions
prop_SEC=rbind(sum_SEC_D,sum_SEC_R)
str(prop_SEC)
prop_SEC$Molecule=as.factor(prop_SEC$Molecule)
prop_SEC$Molecule=factor(prop_SEC$Molecule,levels=rev(levels(prop_SEC$Molecule)))
#Create plot
sec_prop=ggplot(data=prop_SEC,aes(x=reorder(Type,prop),y=prop,fill=Molecule))+
  geom_bar(stat="identity",
           position = position_dodge(),
           color="#333333",
           size=1)+
  geom_errorbar(aes(x=Type,ymin=lower,ymax=upper),width=0.4,colour="#333333",alpha=0.8,size=1,
                position = position_dodge(0.9))+
  theme_classic()+
  scale_fill_manual(name = "Read Set", labels = c("DNA", "RNA"),
                    values=c("#455765","#cd9fb2"))+
  #guide=guide_legend(reverse=T))+
  scale_color_manual(values=c("#455765","#cd9fb2"))+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_x_discrete(name="Gene Cluster Class")+
  scale_y_continuous(name="Proportional Abundance")+
  theme(axis.text.x=element_text(vjust=1,angle=45,hjust=1))
sec_prop

#Add common axis grob
library(gridtext)
yleft = richtext_grob("Relative Abundance (TPM)", rot=90)

#Combine into master plot
library(patchwork)
sec_arrange=((overall_box / NRPS_box) | (overall_box_mat / NRPS_box_mat)) / sec_prop
sec_arrange=sec_arrange + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face='bold'))
sec_arrange + plot_layout(widths=c(1,1),heights=c(4,1),guides="collect") & theme(plot.margin=margin(0))

dev.print(png,fil="sec_arrange.png",type="quartz",antialias="default",width=8.5,
          height=10,units="in",res=1500)
