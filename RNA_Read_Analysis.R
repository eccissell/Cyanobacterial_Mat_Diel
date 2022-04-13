library(dplyr)
library(ggplot2)
library(viridis)
library(tidyverse)
library(gdata)
library(reshape2)
library(phyloseq)
library(vegan)
library(conflicted)
library(ggpubr)
conflict_prefer("summarise",'dplyr')
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")

#####DIEL RNA DATA ALIGNED TO MAGS
#We will have to do some merging to get the feature count matrix
#Read in all the counts
rna_diel1=read.xls('/Users/cissell/Desktop/D1_RNA_counts.xlsx')
rna_diel2=read.xls('/Users/cissell/Desktop/D2_RNA_counts.xlsx')
rna_diel3=read.xls('/Users/cissell/Desktop/D3_RNA_counts.xlsx')
rna_diel4=read.xls('/Users/cissell/Desktop/D4_RNA_counts.xlsx')
#Read in the gene id
rna_kegg1=read.xls('/Users/cissell/Desktop/D1_annotations.xlsx')
rna_kegg2=read.xls('/Users/cissell/Desktop/D2_annotations.xlsx')
rna_kegg3=read.xls('/Users/cissell/Desktop/D3_annotations.xlsx')
rna_kegg4=read.xls('/Users/cissell/Desktop/D4_annotations.xlsx')
#Merge counts with gene ids
merged1=merge(rna_diel1,rna_kegg1, by="contig",all.x=T)
merged2=merge(rna_diel2,rna_kegg2, by="contig",all.x=T)
merged3=merge(rna_diel3,rna_kegg3, by="contig",all.x=T)
merged4=merge(rna_diel4,rna_kegg4, by="contig",all.x=T)
#Remove rows without KEGG ID
smerged1=merged1 %>%
  na_if("") %>%
  na.omit
smerged2=merged2 %>%
  na_if("") %>%
  na.omit
smerged3=merged3 %>%
  na_if("") %>%
  na.omit
smerged4=merged4 %>%
  na_if("") %>%
  na.omit
#Drop contig
smerged1=subset(smerged1,select=-c(contig,Length))
smerged2=subset(smerged2,select=-c(contig,Length))
smerged3=subset(smerged3,select=-c(contig,Length))
smerged4=subset(smerged4,select=-c(contig,Length))
#Summarize within KEGG ID
summerged1=smerged1 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
summerged2=smerged2 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
summerged3=smerged3 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
summerged4=smerged4 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
#Merge across gene ids
rna_merge1=merge(summerged1,summerged2,by="kegg_id",all=T)
rna_merge2=merge(summerged3,summerged4,by="kegg_id",all=T)
rna_merge=merge(rna_merge1,rna_merge2,by="kegg_id",all=T)
rna_merge[is.na(rna_merge)]=0

#Make OTU table for phyloseq
#Rownames have to be ID
rownames(rna_merge)<- rna_merge$kegg_id

#Check counts
count_rna <- colSums(rna_merge[, c(2:ncol(rna_merge))])
count_rna

#Create taxonomy table
rna_tax<- rna_merge %>%
  select(kegg_id)
str(rna_tax)

#Make factors not strings
rna_tax<- rna_tax %>%
  mutate_if(is.character, as.factor)
str(rna_tax)

#Rename first column
colnames(rna_tax)[1]<- "kegg_id"
str(rna_tax)
#OTU IDs have to be rownames, just like OTU table
rownames(rna_tax)<- rna_tax$kegg_id

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
rna_merge<- subset(rna_merge, select=-kegg_id)


#Create metadata table
rna_meta=read.xls('/Users/cissell/Desktop/rna_diel_meta.xlsx')
str(rna_meta)
rownames(rna_meta)<- rna_meta$SampleID
rna_meta<- rna_meta %>%
  select(-SampleID)
#Make proper levels
rna_meta$Mat<- factor(rna_meta$Mat, levels = c("R1", "R2", "R3","R4"))
rna_meta$Age<- factor(rna_meta$Age, levels = c("1", "2","3","4", "5"))

#Create matrices before creating phyloseq object
rna_otu_mat<- as.matrix(rna_merge)
rna_tax_mat<- tax_table(as.matrix(rna_tax))

#Transform to phyloseq objects
rna_phylo_OTU<- otu_table(rna_otu_mat, taxa_are_rows = TRUE)
rna_phylo_TAX<- tax_table(rna_tax_mat)
rna_phylo_samples<- sample_data(rna_meta)

#And unite them into one object
rna_phylo_object<- phyloseq(rna_phylo_OTU, rna_phylo_TAX, rna_phylo_samples)
#Check and see if it makes sense
sample_sums(rna_phylo_object)
sample_names(rna_phylo_object)
rank_names(rna_phylo_object)
sample_variables(rna_phylo_object)
otu_table(rna_phylo_object)[1:3, 1:2]
taxa_names(rna_phylo_object)[1:5]
tax_table(rna_phylo_object)
taxa_names(rna_phylo_object)

library(microbiome)
rna_clr=microbiome::transform(rna_phylo_object, transform='clr',target="sample")
#NMDS on raw dataframe
#MAke dataframe
rna_nmds_df=psmelt(rna_clr)
rna_nmds_df=rna_nmds_df%>%
  select(-OTU)%>%
  spread(rna_nmds_df,key=kegg_id,value=Abundance)

dat=rna_nmds_df[,4:7367]
dat=as.data.frame(lapply(dat,unlist))
str(dat)
grp=as.data.frame(rna_nmds_df[,2:3],header=FALSE)
rna_nmds=metaMDS(dat,distance="euclidean")
rna_nmds
#Stressplot
stressplot(rna_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)
dev.print(png,fil="stress_rna_full.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
#Scree plot of stress
library(goeveg)
dimcheckMDS(
  dat,
  distance = "euclidean",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

kegg.fit=envfit(rna_nmds,dat,permutations=999)
plot(rna_nmds,type="t")
data.score=as.data.frame(scores(rna_nmds))
data.score$site=rownames(data.score)
data.score$grp=grp

species.scores=envfit(rna_nmds,dat)
species.scores2=data.frame((species.scores$vectors)$arrows,(species.scores$vectors)$r,(species.scores$vectors)$pvals)
species.scores2=species.scores2[,-c(3)]
names(species.scores2)[names(species.scores2)=="X.species.scores.vectors..pvals"]="pval"
head(species.scores2)
sig.sp.sc=subset(species.scores2,pval<=0.05)
sig.sp.sc=rownames_to_column(sig.sp.sc, "KEGG")


#Make hull vectors and hull data frame for plotting polygons
age.1 = data.score[(data.score$grp)$Age == "1",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.2 = data.score[(data.score$grp)$Age == "2",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.3=data.score[(data.score$grp)$Age == "3",][chull(data.score[(data.score$grp)$Age== "3",c("NMDS1","NMDS2")]),]
age.4=data.score[(data.score$grp)$Age == "4",][chull(data.score[(data.score$grp)$Age== "4",c("NMDS1","NMDS2")]),]
age.5=data.score[(data.score$grp)$Age == "5",][chull(data.score[(data.score$grp)$Age=="5",c("NMDS1","NMDS2")]),]
hull.data=rbind(age.1,age.2,age.3,age.4,age.5)
str(hull.data)


mat.1 = data.score[(data.score$grp)$Mat == "R1",][chull(data.score[(data.score$grp)$Mat=="R1",c("NMDS1","NMDS2")]),]
mat.2 = data.score[(data.score$grp)$Mat == "R2",][chull(data.score[(data.score$grp)$Mat=="R2",c("NMDS1","NMDS2")]),]
mat.3=data.score[(data.score$grp)$Mat == "R3",][chull(data.score[(data.score$grp)$Mat== "R3",c("NMDS1","NMDS2")]),]
mat.4=data.score[(data.score$grp)$Mat == "R4",][chull(data.score[(data.score$grp)$Mat== "R4",c("NMDS1","NMDS2")]),]
mat.data=rbind(mat.1,mat.2,mat.3,mat.4)
str(mat.data)
#Plot
nmds_rna=ggplot() +
  #geom_polygon(data=hull.data,
  #             aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
  #             alpha=0.4)+
  geom_polygon(data=mat.data,
               aes(x=NMDS1,y=NMDS2, group=grp$Mat),
               #size=1,
               #linetype=1,
               fill="#333333",
               #colour="#333333",
               alpha=0.6)+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2"),name="Time") +
#  scale_fill_manual(values=c("#90a295","#fde992","#455765","#add8e6","#cd9fb2","#fde992")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=100,y=-95, label="Stress = 0.09", color="black",size=5)+
  #annotate(geom="text",x=0,y=-1, label="3", color="#333333",size=5)+
  #annotate(geom="text",x=0.78,y=-0.42, label="1", color="#333333",size=5)+
  #annotate(geom="text",x=-0.8,y=0.2, label="4", color="#333333",size=5)+
  #annotate(geom="text",x=0.3,y=1.25, label="2", color="#333333",size=5)+
  scale_x_continuous(name="NMDS1")+
  scale_y_continuous(name="NMDS2")
nmds_rna

dev.print(png,fil="nmds_rna_full.png",type="quartz",antialias="default",width=6.5,
          height=4.5,units="in",res=1300)

#Create distance matrix from clr
clr_dist_matrix <- phyloseq::distance(rna_clr, method = "euclidean")
distance(phylo_clr, method = "euclidean")

#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Age)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.8047

#Mat
dispr2 <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.9847
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(clr_dist_matrix ~sample_data(rna_clr)$Age+sample_data(rna_clr)$Mat,permutations=999, method="euclidean", by="margin")
permanova

#Pairwise PERMANOVA
#Time
pair_perm=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Age, sim.method="euclidean")
pair_perm #none
#Mat
pair_perm2=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Mat, sim.method="euclidean")
pair_perm2 #2v4 and 3v4


library(MicEco)
#Time
venn_age=ps_venn(rna_phylo_object,group='Age',type='counts',relative=FALSE,fraction=0,fill=c("#90a295","#fde992","#455765","#add8e6","#cd9fb2","#fde992"))
#Individual
venn_mat=ps_venn(rna_phylo_object,group='Mat',type='percent',relative=FALSE,fraction=0,fill=c("#90a295","#455765","#cd9fb2","#fde992"))
#Combine
venns=ggarrange(venn_age,venn_mat,
                ncol=2,nrow=1,align="hv",
                labels=c("a","b"))
venns
annotate_figure(venns,
                left=text_grob("Phylum", color="Black", rot=90,size=15,vjust=(2.3),hjust=0.1))
dev.print(png,fil="venns.png",type="quartz",antialias="default",width=6.5,
          height=4.5,units="in",res=1300)


##################################
#NOW DO IT WITH A SUBSET OF ENERGY FUNCTIONS ONLY, NOT TOTAL PROFILE
##################################
#####DIEL RNA DATA ALIGNED TO MAGS
#We will have to do some merging to get the feature count matrix
#Read in all the counts
rna_diel1=read.xls('/Users/cissell/Desktop/D1_RNA_counts.xlsx')
rna_diel2=read.xls('/Users/cissell/Desktop/D2_RNA_counts.xlsx')
rna_diel3=read.xls('/Users/cissell/Desktop/D3_RNA_counts.xlsx')
rna_diel4=read.xls('/Users/cissell/Desktop/D4_RNA_counts.xlsx')
#Read in the gene id
rna_kegg1=read.xls('/Users/cissell/Desktop/D1_annotations.xlsx')
rna_kegg2=read.xls('/Users/cissell/Desktop/D2_annotations.xlsx')
rna_kegg3=read.xls('/Users/cissell/Desktop/D3_annotations.xlsx')
rna_kegg4=read.xls('/Users/cissell/Desktop/D4_annotations.xlsx')
#Read in energy subsets
energy1=read.xls('/Users/cissell/Desktop/D1_metabolism_summary.xlsx')
energy2=read.xls('/Users/cissell/Desktop/D2_metabolism_summary.xlsx')
energy3=read.xls('/Users/cissell/Desktop/D3_metabolism_summary.xlsx')
energy4=read.xls('/Users/cissell/Desktop/D4_metabolism_summary.xlsx')
#Merge counts with gene ids
merged1=merge(rna_diel1,rna_kegg1, by="contig",all.x=T)
merged2=merge(rna_diel2,rna_kegg2, by="contig",all.x=T)
merged3=merge(rna_diel3,rna_kegg3, by="contig",all.x=T)
merged4=merge(rna_diel4,rna_kegg4, by="contig",all.x=T)
#Remove rows without KEGG ID
smerged1=merged1 %>%
  na_if("") %>%
  na.omit
smerged2=merged2 %>%
  na_if("") %>%
  na.omit
smerged3=merged3 %>%
  na_if("") %>%
  na.omit
smerged4=merged4 %>%
  na_if("") %>%
  na.omit
#Drop contig
smerged1=subset(smerged1,select=-c(contig,Length))
smerged2=subset(smerged2,select=-c(contig,Length))
smerged3=subset(smerged3,select=-c(contig,Length))
smerged4=subset(smerged4,select=-c(contig,Length))
#Summarize within KEGG ID
summerged1=smerged1[1:6] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))
summerged2=smerged2[1:6] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))
summerged3=smerged3[1:6] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))
summerged4=smerged4[1:6] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))

#Merge to subset
enmerged1=merge(energy1,summerged1, by="kegg_id")
enmerged2=merge(energy2,summerged2, by="kegg_id")
enmerged3=merge(energy3,summerged3, by="kegg_id")
enmerged4=merge(energy4,summerged4, by="kegg_id")
#Summarize within KEGG ID
emerged1=enmerged1[c(1,4,5,6,7,8)] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))
emerged2=enmerged2[c(1,4,5,6,7,8)] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))
emerged3=enmerged3[c(1,4,5,6,7,8)] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))
emerged4=enmerged4[c(1,4,5,6,7,8)] %>%
  group_by(kegg_id) %>%
  summarise(across(everything(), sum))
#Merge across gene ids
rna_merge1=merge(y=emerged1,x=emerged2,by="kegg_id",all=T)
rna_merge2=merge(emerged3,emerged4,by="kegg_id",all=T)
rna_merge=merge(y=rna_merge1,x=rna_merge2,by="kegg_id",all=T)

#Replace NA with 0
rna_merge[is.na(rna_merge)] <- 0

#Subset by lower decile
#Create deciles
rna_merge$qR1_1=ntile(rna_merge$R1_1,10)
rna_merge$qR1_2=ntile(rna_merge$R1_2,10)
rna_merge$qR1_3=ntile(rna_merge$R1_3,10)
rna_merge$qR1_4=ntile(rna_merge$R1_4,10)
rna_merge$qR1_5=ntile(rna_merge$R1_5,10)

rna_merge$qR2_1=ntile(rna_merge$R2_1,10)
rna_merge$qR2_2=ntile(rna_merge$R2_2,10)
rna_merge$qR2_3=ntile(rna_merge$R2_3,10)
rna_merge$qR2_4=ntile(rna_merge$R2_4,10)
rna_merge$qR2_5=ntile(rna_merge$R2_5,10)

rna_merge$qR3_1=ntile(rna_merge$R3_1,10)
rna_merge$qR3_2=ntile(rna_merge$R3_2,10)
rna_merge$qR3_3=ntile(rna_merge$R3_3,10)
rna_merge$qR3_4=ntile(rna_merge$R3_4,10)
rna_merge$qR3_5=ntile(rna_merge$R3_5,10)

rna_merge$qR4_1=ntile(rna_merge$R4_1,10)
rna_merge$qR4_2=ntile(rna_merge$R4_2,10)
rna_merge$qR4_3=ntile(rna_merge$R4_3,10)
rna_merge$qR4_4=ntile(rna_merge$R4_4,10)
rna_merge$qR4_5=ntile(rna_merge$R4_5,10)

#Create NA if lower decile
rna_merge$R1_1=with(rna_merge, ifelse(qR1_1 != "1", R1_1, NA))
rna_merge$R1_2=with(rna_merge, ifelse(qR1_2 != "1", R1_2, NA))
rna_merge$R1_3=with(rna_merge, ifelse(qR1_3 != "1", R1_3, NA))
rna_merge$R1_4=with(rna_merge, ifelse(qR1_4 != "1", R1_4, NA))
rna_merge$R1_5=with(rna_merge, ifelse(qR1_5 != "1", R1_5, NA))

rna_merge$R2_1=with(rna_merge, ifelse(qR2_1 != "1", R2_1, NA))
rna_merge$R2_2=with(rna_merge, ifelse(qR2_2 != "1", R2_2, NA))
rna_merge$R2_3=with(rna_merge, ifelse(qR2_3 != "1", R2_3, NA))
rna_merge$R2_4=with(rna_merge, ifelse(qR2_4 != "1", R2_4, NA))
rna_merge$R2_5=with(rna_merge, ifelse(qR2_5 != "1", R2_5, NA))

rna_merge$R3_1=with(rna_merge, ifelse(qR3_1 != "1", R3_1, NA))
rna_merge$R3_2=with(rna_merge, ifelse(qR3_2 != "1", R3_2, NA))
rna_merge$R3_3=with(rna_merge, ifelse(qR3_3 != "1", R3_3, NA))
rna_merge$R3_4=with(rna_merge, ifelse(qR3_4 != "1", R3_4, NA))
rna_merge$R3_5=with(rna_merge, ifelse(qR3_5 != "1", R3_5, NA))

rna_merge$R4_1=with(rna_merge, ifelse(qR4_1 != "1", R4_1, NA))
rna_merge$R4_2=with(rna_merge, ifelse(qR4_2 != "1", R4_2, NA))
rna_merge$R4_3=with(rna_merge, ifelse(qR4_3 != "1", R4_3, NA))
rna_merge$R4_4=with(rna_merge, ifelse(qR4_4 != "1", R4_4, NA))
rna_merge$R4_5=with(rna_merge, ifelse(qR4_5 != "1", R4_5, NA))

#remove deciles
rna_merge=rna_merge[1:21]

#make zero again
rna_merge[is.na(rna_merge)] <- 0

#Make OTU table for phyloseq
#Rownames have to be ID
rownames(rna_merge)<- rna_merge$kegg_id

#Check counts
count_rna <- colSums(rna_merge[, c(2:ncol(rna_merge))])
count_rna

#Create taxonomy table
rna_tax<- rna_merge %>%
  select(kegg_id)
str(rna_tax)

#Make factors not strings
rna_tax<- rna_tax %>%
  mutate_if(is.character, as.factor)
str(rna_tax)

#Rename first column
colnames(rna_tax)[1]<- "kegg_id"
str(rna_tax)
#OTU IDs have to be rownames, just like OTU table
rownames(rna_tax)<- rna_tax$kegg_id

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
rna_merge<- subset(rna_merge, select=-kegg_id)


#Create metadata table
rna_meta=read.xls('/Users/cissell/Desktop/rna_diel_meta.xlsx')
str(rna_meta)
rownames(rna_meta)<- rna_meta$SampleID
rna_meta<- rna_meta %>%
  select(-SampleID)
#Make proper levels
rna_meta$Mat<- factor(rna_meta$Mat, levels = c("R1", "R2", "R3","R4"))
rna_meta$Age<- factor(rna_meta$Age, levels = c("1", "2","3","4", "5"))

#Create matrices before creating phyloseq object
rna_otu_mat<- as.matrix(rna_merge)
rna_tax_mat<- tax_table(as.matrix(rna_tax))

#Transform to phyloseq objects
rna_phylo_OTU<- otu_table(rna_otu_mat, taxa_are_rows = TRUE)
rna_phylo_TAX<- tax_table(rna_tax_mat)
rna_phylo_samples<- sample_data(rna_meta)

#And unite them into one object
rna_phylo_object<- phyloseq(rna_phylo_OTU, rna_phylo_TAX, rna_phylo_samples)
#Check and see if it makes sense
sample_sums(rna_phylo_object)
sample_names(rna_phylo_object)
rank_names(rna_phylo_object)
sample_variables(rna_phylo_object)
otu_table(rna_phylo_object)[1:3, 1:2]
taxa_names(rna_phylo_object)[1:5]
tax_table(rna_phylo_object)
taxa_names(rna_phylo_object)

library(microbiome)
rna_clr=microbiome::transform(rna_phylo_object, transform='clr',target="sample")

#NMDS on raw dataframe
#MAke dataframe
rna_nmds_df=psmelt(rna_clr)
rna_nmds_df=rna_nmds_df%>%
  select(-OTU)%>%
  spread(rna_nmds_df,key=kegg_id,value=Abundance)

dat=rna_nmds_df[,4:272]
dat=as.data.frame(lapply(dat,unlist))
str(dat)
grp=as.data.frame(rna_nmds_df[,2:3],header=FALSE)
rna_nmds=metaMDS(dat,distance="euclidean")
rna_nmds
#Stressplot
stressplot(rna_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)
dev.print(png,fil="stress_rna_subset.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
#Scree plot of stress
library(goeveg)
dimcheckMDS(
  dat,
  distance = "euclidean",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

kegg.fit=envfit(rna_nmds,dat,permutations=999)
plot(rna_nmds,type="t")
data.score=as.data.frame(scores(rna_nmds))
data.score$site=rownames(data.score)
data.score$grp=grp

species.scores=envfit(rna_nmds,dat,permutations=999)
species.scores2=as.data.frame(scores(species.scores,display="vectors"))
species.scores2=cbind(species.scores2, kegg_id=rownames(species.scores2))
species.scores2=cbind(species.scores2,pval=species.scores$vectors$pvals)
head(species.scores2)
sig.sp.sc=subset(species.scores2,pval<=0.05)

#Try to merge with original data frame to get gene functions
sig_drivers_merge=merge(energy1[c(1,3)],sig.sp.sc, by="kegg_id")
sig_drivers_merge=distinct(sig_drivers_merge,kegg_id, .keep_all=T)
conflict_prefer("arrange", "dplyr")
sig_drivers_merge=arrange(sig_drivers_merge,module)

#Create deciles of importance to unnderstand biggest drivers
sig_drivers_merge$abnm1=abs(sig_drivers_merge$NMDS1)
sig_drivers_merge$abnm2=abs(sig_drivers_merge$NMDS2)
sig_drivers_merge$q1=ntile(sig_drivers_merge$abnm1,4)
sig_drivers_merge$q2=ntile(sig_drivers_merge$abnm2,4)
sig_drivers_subset=sig_drivers_merge[(sig_drivers_merge$q1==4) | (sig_drivers_merge$q2==4),]
sig_drivers_subset$NMDS1=sig_drivers_subset$NMDS1*8
sig_drivers_subset$NMDS2=sig_drivers_subset$NMDS2*8

#Add id to rows for subsetting ease
sig_drivers_subset$uid=1:nrow(sig_drivers_subset)
#Subset for important drivers by type
sig_subset_arnon=sig_drivers_subset[c(52:56),]
sig_subset_denitrification=sig_drivers_subset[c(10:11),]
sig_subset_nitrate_red=sig_drivers_subset[c(2:3,14),]
sig_subset_sulf_red=sig_drivers_subset[c(4:6),]
sig_subset_photo1=sig_drivers_subset[c(48:49),]

#Make hull vectors and hull data frame for plotting polygons
age.1 = data.score[(data.score$grp)$Age == "1",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.2 = data.score[(data.score$grp)$Age == "2",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.3=data.score[(data.score$grp)$Age == "3",][chull(data.score[(data.score$grp)$Age== "3",c("NMDS1","NMDS2")]),]
age.4=data.score[(data.score$grp)$Age == "4",][chull(data.score[(data.score$grp)$Age== "4",c("NMDS1","NMDS2")]),]
age.5=data.score[(data.score$grp)$Age == "5",][chull(data.score[(data.score$grp)$Age=="5",c("NMDS1","NMDS2")]),]
hull.data=rbind(age.1,age.2,age.3,age.4,age.5)
str(hull.data)


mat.1 = data.score[(data.score$grp)$Mat == "R1",][chull(data.score[(data.score$grp)$Mat=="R1",c("NMDS1","NMDS2")]),]
mat.2 = data.score[(data.score$grp)$Mat == "R2",][chull(data.score[(data.score$grp)$Mat=="R2",c("NMDS1","NMDS2")]),]
mat.3=data.score[(data.score$grp)$Mat == "R3",][chull(data.score[(data.score$grp)$Mat== "R3",c("NMDS1","NMDS2")]),]
mat.4=data.score[(data.score$grp)$Mat == "R4",][chull(data.score[(data.score$grp)$Mat== "R4",c("NMDS1","NMDS2")]),]
mat.data=rbind(mat.1,mat.2,mat.3,mat.4)
str(mat.data)
#Plot
nmds_rna=ggplot() +
  geom_polygon(data=hull.data,
               aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
               alpha=0.4)+
  geom_polygon(data=mat.data,
               aes(x=NMDS1,y=NMDS2, group=grp$Mat),
               size=1,
               linetype=3,
               fill=NA,
               colour="#333333",
               alpha=0.6)+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2"),name="Time") +
    scale_fill_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2","#fde992")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=-20,y=-12, label="Stress = 0.120", color="black",size=5)+
  #annotate(geom="text",x=0,y=-1, label="3", color="#333333",size=5)+
  #annotate(geom="text",x=0.78,y=-0.42, label="1", color="#333333",size=5)+
  #annotate(geom="text",x=-0.8,y=0.2, label="4", color="#333333",size=5)+
  #annotate(geom="text",x=0.3,y=1.25, label="2", color="#333333",size=5)+
  scale_x_continuous(name="")+
  scale_y_continuous(name="NMDS2")+
  theme(legend.position = "none")
nmds_rna
dev.print(png,fil="nmds_rna_subset.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

nmds_rna2=ggplot() +
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(values=c("#333333","#333333","#333333","#333333","#333333"),name="Time") +
  geom_segment(data=sig_subset_photo1,aes(x=0,y=0,xend=1.5*NMDS1,yend=1.5*NMDS2),size=1.2,alpha=0.8,color="#455765")+
  geom_text(aes(x=10,y=3),label="Photosystem I",color="#455765",size=4,fontface="bold")+
  geom_segment(data=sig_subset_denitrification,aes(x=0,y=0,xend=1.5*NMDS1,yend=1.5*NMDS2),size=1.2,alpha=0.8,color="#90a295")+
  geom_text(aes(x=10,y=-8.5),label="Denitrification",color="#90a295",size=4,fontface="bold")+
  geom_segment(data=sig_subset_nitrate_red,aes(x=0,y=0,xend=1.5*NMDS1,yend=1.5*NMDS2),size=1.2,alpha=0.8,color="#cd9fb2")+
  geom_text(aes(x=-12.5,y=-7),label="Nitrate Reduction",color="#cd9fb2",size=4,fontface="bold")+
  geom_segment(data=sig_subset_sulf_red,aes(x=0,y=0,xend=1.5*NMDS1,yend=1.5*NMDS2),size=1.2,alpha=0.8,color="#add8e6")+
  geom_text(aes(x=-10,y=9),label="Assim. Sulf. Red.",color="#add8e6",size=4,fontface="bold")+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=-20,y=-12, label="Stress = 0.120", color="black",size=5)+
  scale_x_continuous(name="NMDS1")+
  scale_y_continuous(name="NMDS2")+
  theme(legend.position = "none")+
  scale_shape_manual(values=c(19,17,15,18),name="Mat")
nmds_rna2
dev.print(png,fil="nmds_rna_subset_text.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

library(patchwork)
rna_nmds1=nmds_rna/nmds_rna2
rna_nmds1
rna_nmds_arrange=nmds_rna + nmds_rna3 + nmds_rna2 + nmds_rna4
rna_nmds_arrange + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face='bold'))


dev.print(png,fil="rna_nmds_stack.png",type="quartz",antialias="default",width=6.5,height=8.5,units="in",res=1300)

#Create distance matrix from clr
clr_dist_matrix <- phyloseq::distance(rna_clr, method = "euclidean")
#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Age)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.8047

#Mat
dispr2 <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.9844
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(clr_dist_matrix ~sample_data(rna_clr)$Age+sample_data(rna_clr)$Mat,permutations=999, method="euclidean", by="margin")
permanova

#Pairwise PERMANOVA
#Time
pair_perm=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Age, sim.method="euclidean")
pair_perm #none
#Mat
pair_perm2=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Mat, sim.method="euclidean")
pair_perm2 #2v4 and 3v4


library(MicEco)
#Time
venn_age=ps_venn(rna_phylo_object,group='Age',type='counts',relative=FALSE,fraction=0,fill=c("#90a295","#fde992","#455765","#add8e6","#cd9fb2","#fde992"))
#Individual
venn_mat=ps_venn(rna_phylo_object,group='Mat',type='counts',relative=FALSE,fraction=0,fill=c("#90a295","#455765","#cd9fb2","#fde992"))
#Combine
venns=ggarrange(venn_age,venn_mat,
                ncol=2,nrow=1,align="hv",
                labels=c("a","b"))
venns
dev.print(png,fil="venns_RNA_subset.png",type="quartz",antialias="default",width=6.5,
          height=4.5,units="in",res=1300)


##########
##########
##########
#NOW DO SHARED SUBSET OF METABOLIC EXPRESSION FOR TEMPORAL DRIVERS
#Merge across gene ids
rna_merge1=merge(x=emerged1,y=emerged2,by="kegg_id")
rna_merge2=merge(emerged3,emerged4,by="kegg_id")
rna_merge=merge(rna_merge1,rna_merge2,by="kegg_id")
#Replace NA with 0
rna_merge[is.na(rna_merge)] <- 0

#Subset by lower decile
#Create deciles
rna_merge$qR1_1=ntile(rna_merge$R1_1,10)
rna_merge$qR1_2=ntile(rna_merge$R1_2,10)
rna_merge$qR1_3=ntile(rna_merge$R1_3,10)
rna_merge$qR1_4=ntile(rna_merge$R1_4,10)
rna_merge$qR1_5=ntile(rna_merge$R1_5,10)

rna_merge$qR2_1=ntile(rna_merge$R2_1,10)
rna_merge$qR2_2=ntile(rna_merge$R2_2,10)
rna_merge$qR2_3=ntile(rna_merge$R2_3,10)
rna_merge$qR2_4=ntile(rna_merge$R2_4,10)
rna_merge$qR2_5=ntile(rna_merge$R2_5,10)

rna_merge$qR3_1=ntile(rna_merge$R3_1,10)
rna_merge$qR3_2=ntile(rna_merge$R3_2,10)
rna_merge$qR3_3=ntile(rna_merge$R3_3,10)
rna_merge$qR3_4=ntile(rna_merge$R3_4,10)
rna_merge$qR3_5=ntile(rna_merge$R3_5,10)

rna_merge$qR4_1=ntile(rna_merge$R4_1,10)
rna_merge$qR4_2=ntile(rna_merge$R4_2,10)
rna_merge$qR4_3=ntile(rna_merge$R4_3,10)
rna_merge$qR4_4=ntile(rna_merge$R4_4,10)
rna_merge$qR4_5=ntile(rna_merge$R4_5,10)

#Create NA if lower decile
rna_merge$R1_1=with(rna_merge, ifelse(qR1_1 != "1", R1_1, NA))
rna_merge$R1_2=with(rna_merge, ifelse(qR1_2 != "1", R1_2, NA))
rna_merge$R1_3=with(rna_merge, ifelse(qR1_3 != "1", R1_3, NA))
rna_merge$R1_4=with(rna_merge, ifelse(qR1_4 != "1", R1_4, NA))
rna_merge$R1_5=with(rna_merge, ifelse(qR1_5 != "1", R1_5, NA))

rna_merge$R2_1=with(rna_merge, ifelse(qR2_1 != "1", R2_1, NA))
rna_merge$R2_2=with(rna_merge, ifelse(qR2_2 != "1", R2_2, NA))
rna_merge$R2_3=with(rna_merge, ifelse(qR2_3 != "1", R2_3, NA))
rna_merge$R2_4=with(rna_merge, ifelse(qR2_4 != "1", R2_4, NA))
rna_merge$R2_5=with(rna_merge, ifelse(qR2_5 != "1", R2_5, NA))

rna_merge$R3_1=with(rna_merge, ifelse(qR3_1 != "1", R3_1, NA))
rna_merge$R3_2=with(rna_merge, ifelse(qR3_2 != "1", R3_2, NA))
rna_merge$R3_3=with(rna_merge, ifelse(qR3_3 != "1", R3_3, NA))
rna_merge$R3_4=with(rna_merge, ifelse(qR3_4 != "1", R3_4, NA))
rna_merge$R3_5=with(rna_merge, ifelse(qR3_5 != "1", R3_5, NA))

rna_merge$R4_1=with(rna_merge, ifelse(qR4_1 != "1", R4_1, NA))
rna_merge$R4_2=with(rna_merge, ifelse(qR4_2 != "1", R4_2, NA))
rna_merge$R4_3=with(rna_merge, ifelse(qR4_3 != "1", R4_3, NA))
rna_merge$R4_4=with(rna_merge, ifelse(qR4_4 != "1", R4_4, NA))
rna_merge$R4_5=with(rna_merge, ifelse(qR4_5 != "1", R4_5, NA))


#remove deciles
rna_merge=rna_merge[1:21]

#make zero again
rna_merge[is.na(rna_merge)] <- 0

#Make OTU table for phyloseq
#Rownames have to be ID
rownames(rna_merge)<- rna_merge$kegg_id

#Check counts
count_rna <- colSums(rna_merge[, c(2:ncol(rna_merge))])
count_rna

#Create taxonomy table
rna_tax<- rna_merge %>%
  select(kegg_id)
str(rna_tax)

#Make factors not strings
rna_tax<- rna_tax %>%
  mutate_if(is.character, as.factor)
str(rna_tax)

#Rename first column
colnames(rna_tax)[1]<- "kegg_id"
str(rna_tax)
#OTU IDs have to be rownames, just like OTU table
rownames(rna_tax)<- rna_tax$kegg_id

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
rna_merge<- subset(rna_merge, select=-kegg_id)


#Create metadata table
rna_meta=read.xls('/Users/cissell/Desktop/rna_diel_meta.xlsx')
str(rna_meta)
rownames(rna_meta)<- rna_meta$SampleID
rna_meta<- rna_meta %>%
  select(-SampleID)
#Make proper levels
rna_meta$Mat<- factor(rna_meta$Mat, levels = c("R1", "R2", "R3","R4"))
rna_meta$Age<- factor(rna_meta$Age, levels = c("1", "2","3","4", "5"))

#Create matrices before creating phyloseq object
rna_otu_mat<- as.matrix(rna_merge)
rna_tax_mat<- tax_table(as.matrix(rna_tax))

#Transform to phyloseq objects
rna_phylo_OTU<- otu_table(rna_otu_mat, taxa_are_rows = TRUE)
rna_phylo_TAX<- tax_table(rna_tax_mat)
rna_phylo_samples<- sample_data(rna_meta)

#And unite them into one object
rna_phylo_object<- phyloseq(rna_phylo_OTU, rna_phylo_TAX, rna_phylo_samples)
#Check and see if it makes sense
sample_sums(rna_phylo_object)
sample_names(rna_phylo_object)
rank_names(rna_phylo_object)
sample_variables(rna_phylo_object)
otu_table(rna_phylo_object)[1:3, 1:2]
taxa_names(rna_phylo_object)[1:5]
tax_table(rna_phylo_object)
taxa_names(rna_phylo_object)

library(microbiome)
rna_clr=microbiome::transform(rna_phylo_object, transform='clr',target="sample")

#NMDS on raw dataframe
#MAke dataframe
rna_nmds_df=psmelt(rna_clr)
rna_nmds_df=rna_nmds_df%>%
  select(-OTU)%>%
  spread(rna_nmds_df,key=kegg_id,value=Abundance)

dat=rna_nmds_df[,4:210]
dat=as.data.frame(lapply(dat,unlist))
str(dat)
grp=as.data.frame(rna_nmds_df[,2:3],header=FALSE)
rna_nmds=metaMDS(dat,distance="euclidean")
rna_nmds
#Stressplot
stressplot(rna_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)
dev.print(png,fil="stress_rna_subset.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
#Scree plot of stress
library(goeveg)
dimcheckMDS(
  dat,
  distance = "euclidean",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

kegg.fit=envfit(rna_nmds,dat,permutations=999)
plot(rna_nmds,type="t")
data.score=as.data.frame(scores(rna_nmds))
data.score$site=rownames(data.score)
data.score$grp=grp

species.scores=envfit(rna_nmds,dat,permutations=999)
species.scores2=as.data.frame(scores(species.scores,display="vectors"))
species.scores2=cbind(species.scores2, kegg_id=rownames(species.scores2))
species.scores2=cbind(species.scores2,pval=species.scores$vectors$pvals)
head(species.scores2)
sig.sp.sc=subset(species.scores2,pval<=0.05)

#Try to merge with original data frame to get gene functions
sig_drivers_merge=merge(energy1[c(1,3)],sig.sp.sc, by="kegg_id")
sig_drivers_merge=distinct(sig_drivers_merge,kegg_id, .keep_all=T)
conflict_prefer("arrange", "dplyr")
sig_drivers_merge=arrange(sig_drivers_merge,module)
sig_drivers_merge$NMDS1=sig_drivers_merge$NMDS1*10
sig_drivers_merge$NMDS2=sig_drivers_merge$NMDS2*10

#Export csv of drivers
write.csv(sig_drivers_merge,"/Users/cissell/Downloads/share_subset_drivers.csv")
write.csv(sig_drivers_subset,"/Users/cissell/Downloads/subset_drivers.csv")
#Subset for important drivers by type
sig_subset_asim_sulfur=sig_drivers_merge[c(4:7),]
sig_subset_dis_sulfur=sig_drivers_merge[c(37:40),]
sig_subset_nfix=sig_drivers_merge[c(77:79),]
sig_subset_photo2=sig_drivers_merge[c(88:93),]
sig_subset_dis_nitr_red=sig_drivers_merge[c(35:36),]


#Make hull vectors and hull data frame for plotting polygons
age.1 = data.score[(data.score$grp)$Age == "1",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.2 = data.score[(data.score$grp)$Age == "2",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.3=data.score[(data.score$grp)$Age == "3",][chull(data.score[(data.score$grp)$Age== "3",c("NMDS1","NMDS2")]),]
age.4=data.score[(data.score$grp)$Age == "4",][chull(data.score[(data.score$grp)$Age== "4",c("NMDS1","NMDS2")]),]
age.5=data.score[(data.score$grp)$Age == "5",][chull(data.score[(data.score$grp)$Age=="5",c("NMDS1","NMDS2")]),]
hull.data=rbind(age.1,age.2,age.3,age.4,age.5)
str(hull.data)


mat.1 = data.score[(data.score$grp)$Mat == "R1",][chull(data.score[(data.score$grp)$Mat=="R1",c("NMDS1","NMDS2")]),]
mat.2 = data.score[(data.score$grp)$Mat == "R2",][chull(data.score[(data.score$grp)$Mat=="R2",c("NMDS1","NMDS2")]),]
mat.3=data.score[(data.score$grp)$Mat == "R3",][chull(data.score[(data.score$grp)$Mat== "R3",c("NMDS1","NMDS2")]),]
mat.4=data.score[(data.score$grp)$Mat == "R4",][chull(data.score[(data.score$grp)$Mat== "R4",c("NMDS1","NMDS2")]),]
mat.data=rbind(mat.1,mat.2,mat.3,mat.4)
str(mat.data)
#Plot
nmds_rna3=ggplot() +
  geom_polygon(data=hull.data,
               aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
               alpha=0.4)+
  geom_polygon(data=mat.data,
               aes(x=NMDS1,y=NMDS2, group=grp$Mat),
               size=1,
               linetype=3,
               fill=NA,
               colour="#333333",
               alpha=0.6)+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2"),name="Time",labels=c("09:00","15:00","21:00","03:00","09:00")) +
  scale_fill_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2","#fde992"),labels=c("09:00","15:00","21:00","03:00","09:00")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=8,y=-12, label="Stress = 0.177", color="black",size=5)+
  #annotate(geom="text",x=0,y=-1, label="3", color="#333333",size=5)+
  #annotate(geom="text",x=0.78,y=-0.42, label="1", color="#333333",size=5)+
  #annotate(geom="text",x=-0.8,y=0.2, label="4", color="#333333",size=5)+
  #annotate(geom="text",x=0.3,y=1.25, label="2", color="#333333",size=5)+
  scale_x_continuous(name="")+
  scale_y_continuous(name="")
nmds_rna3
dev.print(png,fil="nmds_rna_subset.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

nmds_rna4=ggplot() +
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(values=c("#333333","#333333","#333333","#333333","#333333"),name="Time") +
  geom_segment(data=sig_subset_photo2,aes(x=0,y=0,xend=NMDS1,yend=NMDS2),size=1.2,alpha=0.8,color="#cd9fb2")+
  geom_text(aes(x=8,y=-2),label="Photosystem II",color="#cd9fb2",size=4,fontface="bold")+
  geom_segment(data=sig_subset_dis_nitr_red,aes(x=0,y=0,xend=NMDS1,yend=NMDS2),size=1.2,alpha=0.8,color="#90a295")+
  geom_text(aes(x=2,y=7),label="Dissim. Nitrate Red.",color="#90a295",size=4,fontface="bold")+
  geom_segment(data=sig_subset_nfix,aes(x=0,y=0,xend=NMDS1,yend=NMDS2),size=1.2,alpha=0.8,color="#9c8cdb")+
  geom_text(aes(x=-10,y=1.4),label="Nitrogen Fixation",color="#9c8cdb",size=4,fontface="bold")+
  geom_segment(data=sig_subset_asim_sulfur,aes(x=0,y=0,xend=NMDS1,yend=NMDS2),size=1.2,alpha=0.8,color="#add8e6")+
  geom_text(aes(x=5,y=-6),label="Assim. Sulf. Red.",color="#add8e6",size=4,fontface="bold")+
  geom_segment(data=sig_subset_dis_sulfur,aes(x=0,y=0,xend=NMDS1,yend=NMDS2),size=1.2,alpha=0.8,color="#455765")+
  geom_text(aes(x=-8,y=7),label="Dissim. Sulf. Red.",color="#455765",size=4,fontface="bold")+
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=8,y=-12, label="Stress = 0.177", color="black",size=5)+
  scale_x_continuous(name="NMDS1")+
  scale_y_continuous(name="")+
  theme(legend.position = "none")
nmds_rna4


library(patchwork)
rna_nmds2=nmds_rna3/nmds_rna4
rna_nmds2

rna_nmds_merge=nmds_rna+nmds_rna3+nmds_rna2+nmds_rna4
rna_nmds_merge
rna_nmds_merge + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face='bold'))
dev.print(png,fil="rna_nmds_merge.png",type="quartz",antialias="default",width=10.5,
          height=6.5,units="in",res=1500)

#Create distance matrix from clr
clr_dist_matrix <- phyloseq::distance(rna_clr, method = "euclidean")
#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Age)
dispr
plot(dispr)
permutest(dispr,permutations=20000) #p=0.8047

#Mat
dispr2 <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=20000) #p=0.9844
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(clr_dist_matrix ~sample_data(rna_clr)$Age+sample_data(rna_clr)$Mat,permutations=20000, method="euclidean", by="margin")
permanova

#Pairwise PERMANOVA
#Time
pair_perm=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Age, sim.method="euclidean")
pair_perm #none
#Mat
pair_perm2=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Mat, sim.method="euclidean")
pair_perm2 #2v4 and 3v4
#################
#################
#################
#################
#################
#################
#################
####TPM plots for interpretable visualization of expression
#Merge counts with gene ids
merged1=merge(rna_diel1,rna_kegg1, by="contig",all=T)
merged2=merge(rna_diel2,rna_kegg2, by="contig",all=T)
merged3=merge(rna_diel3,rna_kegg3, by="contig",all=T)
merged4=merge(rna_diel4,rna_kegg4, by="contig",all=T)
#Remove rows without KEGG ID
smerged1=merged1 %>%
  na_if("") %>%
  na.omit
smerged2=merged2 %>%
  na_if("") %>%
  na.omit
smerged3=merged3 %>%
  na_if("") %>%
  na.omit
smerged4=merged4 %>%
  na_if("") %>%
  na.omit
#Drop contig
smerged1=subset(smerged1,select=-c(contig))
smerged2=subset(smerged2,select=-c(contig))
smerged3=subset(smerged3,select=-c(contig))
smerged4=subset(smerged4,select=-c(contig))

#Merge to subset
tpmmerged1=merge(energy1,smerged1, by="kegg_id")
tpmmerged2=merge(energy2,smerged2, by="kegg_id")
tpmmerged3=merge(energy3,smerged3, by="kegg_id")
tpmmerged4=merge(energy4,smerged4, by="kegg_id")

#Add back in SQR
#Subset out by kegg
sqr1=smerged1%>%
  filter(kegg_id=="K17218")
sqr2=smerged2%>%
  filter(kegg_id=="K17218")
sqr3=smerged3%>%
  filter(kegg_id=="K17218")
sqr4=smerged4%>%
  filter(kegg_id=="K17218")
#then combine
tpmmerged1=rbind.fill(tpmmerged1,sqr1)
tpmmerged2=rbind.fill(tpmmerged2,sqr2)
tpmmerged3=rbind.fill(tpmmerged3,sqr3)
tpmmerged4=rbind.fill(tpmmerged4,sqr4)

#Make new df
rpk1=tpmmerged1
rpk2=tpmmerged2
rpk3=tpmmerged3
rpk4=tpmmerged4
#Put length in KB
rpk1$Length=rpk1$Length/1000
rpk2$Length=rpk2$Length/1000
rpk3$Length=rpk3$Length/1000
rpk4$Length=rpk4$Length/1000
#Create RPK for Mat 1
rpk1$R1_1=rpk1$R1_1/rpk1$Length
rpk1$R1_2=rpk1$R1_2/rpk1$Length
rpk1$R1_3=rpk1$R1_3/rpk1$Length
rpk1$R1_4=rpk1$R1_4/rpk1$Length
rpk1$R1_5=rpk1$R1_5/rpk1$Length
#Create RPK for Mat 2
rpk2$R2_1=rpk2$R2_1/rpk2$Length
rpk2$R2_2=rpk2$R2_2/rpk2$Length
rpk2$R2_3=rpk2$R2_3/rpk2$Length
rpk2$R2_4=rpk2$R2_4/rpk2$Length
rpk2$R2_5=rpk2$R2_5/rpk2$Length
#Create RPK for Mat 3
rpk3$R3_1=rpk3$R3_1/rpk3$Length
rpk3$R3_2=rpk3$R3_2/rpk3$Length
rpk3$R3_3=rpk3$R3_3/rpk3$Length
rpk3$R3_4=rpk3$R3_4/rpk3$Length
rpk3$R3_5=rpk3$R3_5/rpk3$Length
#Create RPK for Mat 4
rpk4$R4_1=rpk4$R4_1/rpk4$Length
rpk4$R4_2=rpk4$R4_2/rpk4$Length
rpk4$R4_3=rpk4$R4_3/rpk4$Length
rpk4$R4_4=rpk4$R4_4/rpk4$Length
rpk4$R4_5=rpk4$R4_5/rpk4$Length
#create scaling factors
#Mat 1
sum11=(sum(rpk1$R1_1)/1000000)
sum12=(sum(rpk1$R1_2)/1000000)
sum13=(sum(rpk1$R1_3)/1000000)
sum14=(sum(rpk1$R1_4)/1000000)
sum15=(sum(rpk1$R1_5)/1000000)
#Mat 2
sum21=(sum(rpk2$R2_1)/1000000)
sum22=(sum(rpk2$R2_2)/1000000)
sum23=(sum(rpk2$R2_3)/1000000)
sum24=(sum(rpk2$R2_4)/1000000)
sum25=(sum(rpk2$R2_5)/1000000)
#Mat 3
sum31=(sum(rpk3$R3_1)/1000000)
sum32=(sum(rpk3$R3_2)/1000000)
sum33=(sum(rpk3$R3_3)/1000000)
sum34=(sum(rpk3$R3_4)/1000000)
sum35=(sum(rpk3$R3_5)/1000000)
#Mat 4
sum41=(sum(rpk4$R4_1)/1000000)
sum42=(sum(rpk4$R4_2)/1000000)
sum43=(sum(rpk4$R4_3)/1000000)
sum44=(sum(rpk4$R4_4)/1000000)
sum45=(sum(rpk4$R4_5)/1000000)
#Create TPM for Mat 1
rpk1$R1_1=rpk1$R1_1/sum11
rpk1$R1_2=rpk1$R1_2/sum12
rpk1$R1_3=rpk1$R1_3/sum13
rpk1$R1_4=rpk1$R1_4/sum14
rpk1$R1_5=rpk1$R1_5/sum15
#Create TPM for Mat 2
rpk2$R2_1=rpk2$R2_1/sum21
rpk2$R2_2=rpk2$R2_2/sum22
rpk2$R2_3=rpk2$R2_3/sum23
rpk2$R2_4=rpk2$R2_4/sum24
rpk2$R2_5=rpk2$R2_5/sum25
#Create TPM for Mat 3
rpk3$R3_1=rpk3$R3_1/sum31
rpk3$R3_2=rpk3$R3_2/sum32
rpk3$R3_3=rpk3$R3_3/sum33
rpk3$R3_4=rpk3$R3_4/sum34
rpk3$R3_5=rpk3$R3_5/sum35
#Create TPM for Mat 4
rpk4$R4_1=rpk4$R4_1/sum41
rpk4$R4_2=rpk4$R4_2/sum42
rpk4$R4_3=rpk4$R4_3/sum43
rpk4$R4_4=rpk4$R4_4/sum44
rpk4$R4_5=rpk4$R4_5/sum45

#Check for sense
sum(rpk1$R1_1) #S'all good

#Merge across KEGG ID and Class
tpm1=rpk1 %>%
  group_by(kegg_id,Phylum,Class) %>%
  summarise(across(3:8, sum))%>%
  select(-Length)
tpm2=rpk2 %>%
  group_by(kegg_id,Phylum,Class) %>%
  summarise(across(3:8, sum))%>%
  select(-Length)
tpm3=rpk3 %>%
  group_by(kegg_id,Phylum,Class) %>%
  summarise(across(3:8, sum))%>%
  select(-Length)
tpm4=rpk4 %>%
  group_by(kegg_id,Phylum,Class) %>%
  summarise(across(3:8, sum))%>%
  select(-Length)

#Check again for sense
sum(tpm3$R3_1) #S'all good

#Merge to get descriptions again
tpmmerged1=merge(energy1,tpm1, by="kegg_id",all.y=T)
tpmmerged2=merge(energy2,tpm2, by="kegg_id",all.y=T)
tpmmerged3=merge(energy3,tpm3, by="kegg_id",all.y=T)
tpmmerged4=merge(energy4,tpm4, by="kegg_id",all.y=T)

#Try without getting descriptions again
tpmmerged1=tpm1
tpmmerged2=tpm2
tpmmerged3=tpm3
tpmmerged4=tpm4
  
#Remove redundancy
tpmmerged1=distinct(tpmmerged1,kegg_id,Phylum,Class,.keep_all = T)
tpmmerged2=distinct(tpmmerged2,kegg_id,Phylum,Class,.keep_all = T)
tpmmerged3=distinct(tpmmerged3,kegg_id,Phylum,Class,.keep_all = T)
tpmmerged4=distinct(tpmmerged4,kegg_id,Phylum,Class,.keep_all = T)

#Check again for sense
sum(tpmmerged2$R2_4) #S'all good

#Merge across gene ids
tpm_merge1=merge(y=tpmmerged1,x=tpmmerged2,by=c("kegg_id","Phylum","Class"),all=T)
tpm_merge2=merge(tpmmerged3,tpmmerged4,by=c("kegg_id","Phylum","Class"),all=T)
tpm_merge=merge(y=tpm_merge1,x=tpm_merge2,by=c("kegg_id","Phylum","Class"),all=T)

#Replace NA with 0
tpm_merge[is.na(tpm_merge)] <- 0

#Check for sense
sum(tpm_merge$R3_4)
sum(tpm_merge$R1_1)
#Convert to long
tpm_long <- gather(tpm_merge, Sample, TPM, R3_1:R1_5)

#Add column for time and MatID
conflict_prefer("mutate", "dplyr")
conflict_prefer("startsWith", "gdata")
tpm_long=tpm_long %>%
  mutate(Time=case_when(
    endsWith(Sample,"1") ~ "1",
    endsWith(Sample,"2") ~ "2",
    endsWith(Sample,"3") ~ "3",
    endsWith(Sample,"4") ~ "4",
    endsWith(Sample,"5") ~ "5",
  ))
tpm_long=tpm_long %>%
  mutate(Mat=case_when(
    str_detect(Sample,"R1") ~ "R1",
    str_detect(Sample,"R2") ~ "R2",
    str_detect(Sample,"R3") ~ "R3",
    str_detect(Sample,"R4") ~ "R4",
    str_detect(Sample,"R5") ~ "R5",
  ))%>%
  mutate_at("Class",str_replace,"c__","")%>%
  mutate_at("Class",str_replace,"__F","")%>%
  mutate_at("Phylum",str_replace,"p__","")%>%
  mutate_at("Phylum",str_replace,"_F","")
#Replace blank class with NA, then replace with Phylum
tpm_long=tpm_long%>%
  mutate(Class=na_if(Class,""))%>%
  mutate(Class=ifelse(is.na(Class),Phylum,Class))

#Create stacked area chart for phylum/class across time for all metabolic functions
tpm_long_area_sub=tpm_long[!str_detect(tpm_long$gene_description, "psaA|psaB|psbA|psbD"),]
tpm_stack_df=tpm_long_area_sub%>%
  group_by(Time,Phylum)%>%
  summarise(sumTPM=sum(TPM))%>%
  mutate(percent=sumTPM/sum(sumTPM))
tpm_stack_df=as.data.frame(tpm_stack_df)
tpm_stack_df$Time=as.numeric(tpm_stack_df$Time)
str(tpm_stack_df)
tpm_area_plot=ggplot(data=tpm_stack_df,aes(x=Time,y=percent,fill=Phylum))+
  geom_area(alpha=0.8,size=0.5,colour="white")+
  scale_fill_viridis(discrete=T)+
  theme_classic()+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1))+
  xlab("Sampling time point") +
  ylab("Proportion of metabolic transcripts (TPM)")
tpm_area_plot
dev.print(png,fil="diel_area.png",type="quartz",antialias="default",width=5,
          height=5,units="in",res=1500)
#Subset out important genes (Carbon, Sulfur, Nitrogen, Phosphorous)

#Subset out important genes (Carbon, Sulfur, Nitrogen, Phosphorous)
tpm_long_subs=tpm_long[str_detect(tpm_long$kegg_id,"K02703|K02706|K02705|K02689|K02694|K02690|K02691|K02692|K02274|K02275|K02276|K02277|K17218|K13811|K00958|K00956|K00957|K00860|K00390|K00380|K00381|K00392|K00394|K00395|K11180|K11181|K17222|K17223|K17226|K17227|K17224|K17225|K22622|K02588|K10946|K10944|K10945|K10535|K00371|K00370|K00374|K02567|K02568|K00368|K15864|K04561|K02305|K00376|K00362|K00363|K03385|K15876|K00367|K10534|K00372|K00360|K17877|K00366|K20932|K20933|K20934|K20935|K04758"),]
conflict_prefer("filter", "dplyr")
#Filter out errant classification
tpm_long_subs=tpm_long_subs%>%
  filter(!(kegg_id=="K02705" &
             Class=="Bacteroidia"))%>%
  filter(!(kegg_id=="K02691" &
             Class=="Gammaproteobacteria"))%>%
  filter(!(kegg_id=="K02691" &
             Class=="Bacteroidia"))
tpm_long_subs=tpm_long_subs%>%
  filter(!(kegg_id=="K02706" &
             Class=="Bacteroidia"))%>%
  filter(!(kegg_id=="K02706" &
             Class=="Gammaproteobacteria"))
str(tpm_long_subs)
#Plot
#Set factor order of Kegg ID for plotting
tpm_long_subs$kegg_id=factor(tpm_long_subs$kegg_id,levels=c("K02703","K02706","K02705","K02689","K02694","K02690","K02691","K02692","K02274","K02275","K02276","K02277","K17218","K13811","K00958","K00956","K00957","K00860","K00390","K00380","K00381","K00392","K00394","K00395","K11180","K11181","K17222","K17223","K17226","K17227","K17224","K17225","K22622","K02588","K10946","K10944","K10945","K10535","K00370","K00371","K00374","K02567","K02568","K00368K15864","K04561","K02305","K00376","K00362","K00363","K03385","K15876","K00367","K10534","K00372","K00360","K17877","K00366","K20932","K20933","K20934","K20935","K04758"))

tpm_long_subs=tpm_long_subs[!is.na(tpm_long_subs$kegg_id),]


#Remove 0s from dataframe
tpm_long_subs=tpm_long_subs[tpm_long_subs$TPM !=0,]

#Plot
diel_rna_heat2=ggplot(data=tpm_long_subs,aes(x=Class,
                                          y=kegg_id,
                                          fill=log10(1+TPM)))+
  geom_tile()+
  theme_classic()+
  facet_grid(~Time)+
  scale_fill_viridis_c(option="plasma",
                       na.value="White")+
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1))+
  xlab("")+
  ylab("")+
  labs(fill="Log(TPM)")
diel_rna_heat2
ggsave(diel_rna_heat2,filename="/Users/cissell/Desktop/diel_heats_time.tiff")
#Plot
diel_rna_heat3=ggplot(data=tpm_long_subs,aes(x=Class,
                                             y=kegg_id,
                                             fill=log10(1+TPM)))+
  geom_tile()+
  theme_classic()+
  facet_grid(~Mat)+
  scale_fill_viridis_c(option="plasma",
                       na.value="White")+
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1))+
  theme(axis.text.x=element_blank())+
  xlab("")+
  ylab("")+
  labs(fill="Log(TPM)")
diel_rna_heat3
#Try by maximum expression for each gene_description across time
tpm_long_max=tpm_long_subs%>%
  group_by(kegg_id,Mat,Class)%>%
  mutate(TPM=TPM/max(TPM))

tpm_long_max$TPM[is.na(tpm_long_max$TPM)]=0
tpm_long_max=tpm_long_max[!is.na(tpm_long_max$kegg_id),]

#Now plot maximal
diel_rna_heat_mat=ggplot(data=tpm_long_max,aes(x=Class,
                                            y=kegg_id,
                                            fill=TPM))+
  geom_tile()+
  theme_classic()+
  facet_grid(~Time)+
  scale_fill_viridis_c(na.value="Black")+
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1))+
#  theme(strip.background = element_blank(),
#        strip.text = element_blank())+
  ylab("")+
  labs(fill="Prop. Max TPM")
diel_rna_heat_mat

library(patchwork)
diel_rna_heats_combine=(diel_rna_heat3/diel_rna_heat_mat) + plot_annotation(tag_levels = 'a') &
  theme(plot.tag=element_text(face="bold"))
diel_rna_heats_combine=diel_rna_heats_combine + plot_layout(guides="collect")&theme(legend.position = "right")
diel_rna_heats_combine
ggsave(diel_rna_heats_combine,filename="/Users/cissell/Desktop/diel_heats.tiff")

###############################
###############################
###############################
###############################
###############################
##########################
##########################
##########################
##########################
##########################
##########################
#####COMPARE MAG COMP TO MT AND READ BASED##################
##########################
##########################
#First read in genome taxonomy
GTD_1=read_tsv("/Users/cissell/Desktop/GTD_D1.bac120.summary.tsv")
GTD_2=read_tsv("/Users/cissell/Desktop/GTD_D2.bac120.summary.tsv")
GTD_3=read_tsv("/Users/cissell/Desktop/GTD_D3.bac120.summary.tsv")
GTD_4=read_tsv("/Users/cissell/Desktop/GTD_D4.bac120.summary.tsv")
#Now read in MAG abundances
MAG_1=read_tsv("/Users/cissell/Desktop/D1_genome_props.tsv")
MAG_1=MAG_1[,c(1,8:10)]
MAG_2=read_tsv("/Users/cissell/Desktop/D2_genome_props.tsv")
MAG_2=MAG_2[,c(1,8:10)]
MAG_3=read_tsv("/Users/cissell/Desktop/D3_genome_props.tsv")
MAG_3=MAG_3[,c(1,8:10)]
MAG_4=read_tsv("/Users/cissell/Desktop/D4_genome_props.tsv")
MAG_4=MAG_4[,c(1,8:10)]

#Merge with taxonomy
MAG_merge1=merge(MAG_1,GTD_1, by="Genome",all.x=T)
MAG_merge2=merge(MAG_2,GTD_2, by="Genome",all.x=T)
MAG_merge3=merge(MAG_3,GTD_3, by="Genome",all.x=T)
MAG_merge4=merge(MAG_4,GTD_4, by="Genome",all.x=T)

#Convert wide to long
MAG_long_1=gather(MAG_merge1,Sample, Prop, 2:4)
MAG_long_2=gather(MAG_merge2,Sample, Prop, 2:4)
MAG_long_3=gather(MAG_merge3,Sample, Prop, 2:4)
MAG_long_4=gather(MAG_merge4,Sample, Prop, 2:4)

#Merge into one master frame
MAG_long=rbind(MAG_long_1,MAG_long_2,MAG_long_3,MAG_long_4)

#Add new Time and Mat column
conflict_prefer("mutate", "dplyr")
MAG_long=MAG_long %>%
  mutate(Time=case_when(
    endsWith(Sample,"1") ~ "1",
    endsWith(Sample,"3") ~ "3",
    endsWith(Sample,"5") ~ "5",
  ))
MAG_long=MAG_long %>%
  mutate(Mat=case_when(
    str_detect(Sample,"D1") ~ "1",
    str_detect(Sample,"D2") ~ "2",
    str_detect(Sample,"D3") ~ "3",
    str_detect(Sample,"D4") ~ "4",
    str_detect(Sample,"R5") ~ "R5",
  ))

#Add MAG notation for later plotting of all differnt kinds
MAG_long=cbind(MAG_long,Source='MAG')

#Add in the Read-Based MG data
conflict_prefer("rename", "dplyr")
mg_barplot_df=barplot_df
mg_barplot_df = mg_barplot_df %>%
  rename(
    Prop=Abundance
  )
mg_barplot_df=cbind(mg_barplot_df,Source='MGR')
mg_barplot_df=mg_barplot_df %>%
  rename(
    Time=Age
  )
mg_barplot_df$Mat =gsub('D1','1',mg_barplot_df$Mat)
mg_barplot_df$Mat =gsub('D2','2',mg_barplot_df$Mat)
mg_barplot_df$Mat =gsub('D3','3',mg_barplot_df$Mat)
mg_barplot_df$Mat =gsub('D4','4',mg_barplot_df$Mat)

#Add in total MT
MT_1=merge(rna_diel1,rna_kegg1, by="contig",all.x=T)
MT_2=merge(rna_diel2,rna_kegg2, by="contig",all.x=T)
MT_3=merge(rna_diel3,rna_kegg3, by="contig",all.x=T)
MT_4=merge(rna_diel4,rna_kegg4, by="contig",all.x=T)
#Remove rows without KEGG ID
MT_1=MT_1 %>%
  na_if("") %>%
  na.omit
MT_2=MT_2 %>%
  na_if("") %>%
  na.omit
MT_3=MT_3 %>%
  na_if("") %>%
  na.omit
MT_4=MT_4 %>%
  na_if("") %>%
  na.omit
#Drop contig
MT_1=subset(MT_1,select=-c(contig))
MT_2=subset(MT_2,select=-c(contig))
MT_3=subset(MT_3,select=-c(contig))
MT_4=subset(MT_4,select=-c(contig))
#Make new df
rpk1=MT_1
rpk2=MT_2
rpk3=MT_3
rpk4=MT_4
#Put length in KB
rpk1$Length=rpk1$Length/1000
rpk2$Length=rpk2$Length/1000
rpk3$Length=rpk3$Length/1000
rpk4$Length=rpk4$Length/1000
#Create RPK for Mat 1
rpk1$R1_1=rpk1$R1_1/rpk1$Length
rpk1$R1_2=rpk1$R1_2/rpk1$Length
rpk1$R1_3=rpk1$R1_3/rpk1$Length
rpk1$R1_4=rpk1$R1_4/rpk1$Length
rpk1$R1_5=rpk1$R1_5/rpk1$Length
#Create RPK for Mat 2
rpk2$R2_1=rpk2$R2_1/rpk2$Length
rpk2$R2_2=rpk2$R2_2/rpk2$Length
rpk2$R2_3=rpk2$R2_3/rpk2$Length
rpk2$R2_4=rpk2$R2_4/rpk2$Length
rpk2$R2_5=rpk2$R2_5/rpk2$Length
#Create RPK for Mat 3
rpk3$R3_1=rpk3$R3_1/rpk3$Length
rpk3$R3_2=rpk3$R3_2/rpk3$Length
rpk3$R3_3=rpk3$R3_3/rpk3$Length
rpk3$R3_4=rpk3$R3_4/rpk3$Length
rpk3$R3_5=rpk3$R3_5/rpk3$Length
#Create RPK for Mat 4
rpk4$R4_1=rpk4$R4_1/rpk4$Length
rpk4$R4_2=rpk4$R4_2/rpk4$Length
rpk4$R4_3=rpk4$R4_3/rpk4$Length
rpk4$R4_4=rpk4$R4_4/rpk4$Length
rpk4$R4_5=rpk4$R4_5/rpk4$Length
#create scaling factors
#Mat 1
sum11=(sum(rpk1$R1_1)/1000000)
sum12=(sum(rpk1$R1_2)/1000000)
sum13=(sum(rpk1$R1_3)/1000000)
sum14=(sum(rpk1$R1_4)/1000000)
sum15=(sum(rpk1$R1_5)/1000000)
#Mat 2
sum21=(sum(rpk2$R2_1)/1000000)
sum22=(sum(rpk2$R2_2)/1000000)
sum23=(sum(rpk2$R2_3)/1000000)
sum24=(sum(rpk2$R2_4)/1000000)
sum25=(sum(rpk2$R2_5)/1000000)
#Mat 3
sum31=(sum(rpk3$R3_1)/1000000)
sum32=(sum(rpk3$R3_2)/1000000)
sum33=(sum(rpk3$R3_3)/1000000)
sum34=(sum(rpk3$R3_4)/1000000)
sum35=(sum(rpk3$R3_5)/1000000)
#Mat 4
sum41=(sum(rpk4$R4_1)/1000000)
sum42=(sum(rpk4$R4_2)/1000000)
sum43=(sum(rpk4$R4_3)/1000000)
sum44=(sum(rpk4$R4_4)/1000000)
sum45=(sum(rpk4$R4_5)/1000000)
#Create TPM for Mat 1
rpk1$R1_1=rpk1$R1_1/sum11
rpk1$R1_2=rpk1$R1_2/sum12
rpk1$R1_3=rpk1$R1_3/sum13
rpk1$R1_4=rpk1$R1_4/sum14
rpk1$R1_5=rpk1$R1_5/sum15
#Create TPM for Mat 2
rpk2$R2_1=rpk2$R2_1/sum21
rpk2$R2_2=rpk2$R2_2/sum22
rpk2$R2_3=rpk2$R2_3/sum23
rpk2$R2_4=rpk2$R2_4/sum24
rpk2$R2_5=rpk2$R2_5/sum25
#Create TPM for Mat 3
rpk3$R3_1=rpk3$R3_1/sum31
rpk3$R3_2=rpk3$R3_2/sum32
rpk3$R3_3=rpk3$R3_3/sum33
rpk3$R3_4=rpk3$R3_4/sum34
rpk3$R3_5=rpk3$R3_5/sum35
#Create TPM for Mat 4
rpk4$R4_1=rpk4$R4_1/sum41
rpk4$R4_2=rpk4$R4_2/sum42
rpk4$R4_3=rpk4$R4_3/sum43
rpk4$R4_4=rpk4$R4_4/sum44
rpk4$R4_5=rpk4$R4_5/sum45

tpm1=rpk1 %>%
  group_by(kegg_id,Phylum) %>%
  summarise(across(1:5, sum))
tpm2=rpk2 %>%
  group_by(kegg_id,Phylum) %>%
  summarise(across(1:5, sum))
tpm3=rpk3 %>%
  group_by(kegg_id,Phylum) %>%
  summarise(across(1:5, sum))
tpm4=rpk4 %>%
  group_by(kegg_id,Phylum) %>%
  summarise(across(1:5, sum))

MT_merge1=merge(y=tpm1,x=tpm2,by=c("kegg_id","Phylum"),all.x=T)
MT_merge2=merge(tpm3,tpm4,by=c("kegg_id","Phylum"),all.x=T)
MT_merge=merge(y=MT_merge1,x=MT_merge2,by=c("kegg_id","Phylum"),all.x=T)

#Replace NA with 0
MT_merge[is.na(MT_merge)] <- 0

#Convert to long
MT_long <- gather(MT_merge, Sample, TPM, R3_1:R1_5)

#Add column for time and MatID
conflict_prefer("mutate", "dplyr")
conflict_prefer("startsWith", "gdata")

MT_long=MT_long %>%
  mutate(Time=case_when(
    endsWith(Sample,"1") ~ "1",
    endsWith(Sample,"2") ~ "2",
    endsWith(Sample,"3") ~ "3",
    endsWith(Sample,"4") ~ "4",
    endsWith(Sample,"5") ~ "5",
  ))
MT_long=MT_long %>%
  mutate(Mat=case_when(
    str_detect(Sample,"R1") ~ "1",
    str_detect(Sample,"R2") ~ "2",
    str_detect(Sample,"R3") ~ "3",
    str_detect(Sample,"R4") ~ "4",
    str_detect(Sample,"R5") ~ "5",
  ))

#Summarize to get proportions by phylum
sum_MT=MT_long %>%
  group_by(Phylum, Sample) %>%
  summarise(across(TPM, sum))

sum1=subset(sum_MT,Sample=="R1_1")
sum1$prop=sum1$TPM/sum(sum1$TPM)

sum2=subset(sum_MT,Sample=="R1_2")
sum2$prop=sum2$TPM/sum(sum2$TPM)

sum3=subset(sum_MT,Sample=="R1_3")
sum3$prop=sum3$TPM/sum(sum3$TPM)

sum4=subset(sum_MT,Sample=="R1_4")
sum4$prop=sum4$TPM/sum(sum4$TPM)

sum5=subset(sum_MT,Sample=="R1_5")
sum5$prop=sum5$TPM/sum(sum5$TPM)

sum6=subset(sum_MT,Sample=="R2_1")
sum6$prop=sum6$TPM/sum(sum6$TPM)

sum7=subset(sum_MT,Sample=="R2_2")
sum7$prop=sum7$TPM/sum(sum7$TPM)

sum8=subset(sum_MT,Sample=="R2_3")
sum8$prop=sum8$TPM/sum(sum8$TPM)

sum9=subset(sum_MT,Sample=="R2_4")
sum9$prop=sum9$TPM/sum(sum9$TPM)

sum11=subset(sum_MT,Sample=="R2_5")
sum11$prop=sum11$TPM/sum(sum11$TPM)

sum12=subset(sum_MT,Sample=="R3_1")
sum12$prop=sum12$TPM/sum(sum12$TPM)

sum13=subset(sum_MT,Sample=="R3_2")
sum13$prop=sum13$TPM/sum(sum13$TPM)

sum14=subset(sum_MT,Sample=="R3_3")
sum14$prop=sum14$TPM/sum(sum14$TPM)

sum15=subset(sum_MT,Sample=="R3_4")
sum15$prop=sum15$TPM/sum(sum15$TPM)

sum16=subset(sum_MT,Sample=="R3_5")
sum16$prop=sum16$TPM/sum(sum16$TPM)

sum17=subset(sum_MT,Sample=="R4_1")
sum17$prop=sum17$TPM/sum(sum17$TPM)

sum18=subset(sum_MT,Sample=="R4_2")
sum18$prop=sum18$TPM/sum(sum18$TPM)

sum19=subset(sum_MT,Sample=="R4_3")
sum19$prop=sum19$TPM/sum(sum19$TPM)

sum20=subset(sum_MT,Sample=="R4_4")
sum20$prop=sum20$TPM/sum(sum20$TPM)

sum21=subset(sum_MT,Sample=="R4_5")
sum21$prop=sum21$TPM/sum(sum21$TPM)

MT_prop_merge=rbind(sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum11,sum12,sum13,sum14,sum15,sum16,sum17,sum18,sum19,sum20,sum21)
MT_prop_merge=subset(MT_prop_merge, TPM != "0")
MT_prop_merge$prop=MT_prop_merge$prop*100

#Add in source column and rename prop
MT_prop_merge = MT_prop_merge %>%
  rename(
    Prop=prop
  )
MT_prop_merge=cbind(MT_prop_merge,Source='MT')

#Add new Time and Mat column
MT_prop_merge=MT_prop_merge %>%
  mutate(Time=case_when(
    endsWith(Sample,"1") ~ "1",
    endsWith(Sample,"2") ~ "2",
    endsWith(Sample,"3") ~ "3",
    endsWith(Sample,"4") ~ "4",
    endsWith(Sample,"5") ~ "5",
  ))
MT_prop_merge=MT_prop_merge %>%
  mutate(Mat=case_when(
    str_detect(Sample,"R1") ~ "1",
    str_detect(Sample,"R2") ~ "2",
    str_detect(Sample,"R3") ~ "3",
    str_detect(Sample,"R4") ~ "4",
  ))
#Subset to only have similar times as DNA
MT_prop_sub=MT_prop_merge[!(MT_prop_merge$Time=='2' | MT_prop_merge$Time=='4'),]

##Combine into master df for plotting
abund_master_df=bind_rows(mg_barplot_df,MAG_long,MT_prop_sub)

abund_master_df_sub=abund_master_df[,c(3:7,8)]
abund_master_df_sub$Phylum <- gsub('p__', '', abund_master_df_sub$Phylum)
abund_master_df_sub$Phylum=as.factor(abund_master_df_sub$Phylum)
abund_master_df_sub$Time=as.factor(abund_master_df_sub$Time)
abund_master_df_sub$Mat=as.factor(abund_master_df_sub$Mat)
abund_master_df_sub$Source=as.factor(abund_master_df_sub$Source)
str(abund_master_df_sub)

#Rename some levels for better plotting
abund_master_df_sub$Phylum =gsub('Acidobacteriota','Acidobacteria',abund_master_df_sub$Phylum)
abund_master_df_sub$Phylum =gsub('Bacteroidota','Bacteroidetes',abund_master_df_sub$Phylum)
abund_master_df_sub$Phylum =gsub('Planctomycetota','Planctomycetes',abund_master_df_sub$Phylum)
abund_master_df_sub$Phylum =gsub('Verrucomicrobiota','Verrucomicrobia',abund_master_df_sub$Phylum)

abund_master_df_sub$Time =gsub('1','Time 1 - 09:00',abund_master_df_sub$Time)
abund_master_df_sub$Time =gsub('3','Time 3 - 21:00',abund_master_df_sub$Time)
abund_master_df_sub$Time =gsub('5','Time 5 - 09:00',abund_master_df_sub$Time)
str(abund_master_df_sub)
abund_master_df_sub$Source=forcats::fct_inorder(abund_master_df_sub$Source)
abund_master_df_sub$Mat =gsub('1','Mat 1',abund_master_df_sub$Mat)
abund_master_df_sub$Mat =gsub('2','Mat 2',abund_master_df_sub$Mat)
abund_master_df_sub$Mat =gsub('3','Mat 3',abund_master_df_sub$Mat)
abund_master_df_sub$Mat =gsub('4','Mat 4',abund_master_df_sub$Mat)
#Summarize into mean and SE by mat
master_sum_mat <- abund_master_df_sub %>%
  group_by(Mat, Time, Phylum, Source) %>%
  summarise(Abundsum = sum(Prop))%>%
  group_by(Mat, Phylum, Source) %>%
  summarise(N = n(),
            mAbund = mean(Abundsum),
            seAbund = sd(Abundsum)/sqrt(N))

#Summarize by time mean SE
master_sum_time <- abund_master_df_sub %>%
  group_by(Mat, Time, Phylum, Source) %>%
  dplyr::summarise(Abundsum = sum(Prop))%>%
  group_by(Time, Phylum, Source) %>%
  dplyr::summarise(N = n(),
                   mAbund = mean(Abundsum),
                   seAbund = sd(Abundsum)/sqrt(N))
#Plot
perct_time=ggplot(data=subset(master_sum_time,mAbund>1 & Phylum !="Above Phyla" & Phylum !="<1%"),aes(x=Phylum,y=mAbund,fill=Source)) +
  geom_col(position="dodge",width=0.9) +
  geom_errorbar(aes(ymin=mAbund-seAbund, ymax=mAbund+seAbund),
                position=position_dodge(0.9),
                width=0,
                size=0.8)+
  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2"))+
  labs(y= "Relative abundance [%]",
       fill= "Source")+
  facet_wrap(~Time) +
  theme_classic()+
  ylab("Relative Abundance (%)") +
  xlab("")+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.title.y = element_text(margin =
                                      margin(t = 0, r = 8, b = 0, l = 0)))+
  theme(strip.text.x = element_text(size=13)) +
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
perct_time

perct_mat=ggplot(data=subset(master_sum_mat,mAbund>1 & Phylum !="Above Phyla" & Phylum !="<1%"),aes(x=Phylum,y=mAbund,fill=Source)) +
  geom_col(position="dodge",width=0.9) +
  geom_errorbar(aes(ymin=mAbund-seAbund, ymax=mAbund+seAbund),
                position=position_dodge(0.9),
                width=0,
                size=0.8)+
  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2"))+
  labs(y= "Relative abundance [%]",
       fill= "Source")+
  facet_wrap(~Mat) +
  theme_classic()+
  ylab("Relative Abundance (%)") +
  xlab("Phylum")+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.title.y = element_text(margin =
                                      margin(t = 0, r = 8, b = 0, l = 0))) +
  theme(strip.text.x = element_text(size=13)) +
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
perct_mat

library(patchwork)
master_perct_arrange=perct_time / perct_mat
master_perct_arrange=master_perct_arrange + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face='bold'))
master_perct_arrange + plot_layout(guides="collect")& theme(legend.position = "right")

dev.print(png,fil="master_abund.png",type="quartz",antialias="default",width=8.5,
          height=10,units="in",res=1500)




###################
####################
#Normalize to expression of constituitive genes
#recA; rpoB; proC; secA
####################
###################
CG_1=merge(rna_diel1,rna_kegg1, by="contig",all.x=T)
CG_2=merge(rna_diel2,rna_kegg2, by="contig",all.x=T)
CG_3=merge(rna_diel3,rna_kegg3, by="contig",all.x=T)
CG_4=merge(rna_diel4,rna_kegg4, by="contig",all.x=T)
#Remove rows without KEGG ID
CG_1=CG_1 %>%
  na_if("") %>%
  na.omit
CG_2=CG_2 %>%
  na_if("") %>%
  na.omit
CG_3=CG_3 %>%
  na_if("") %>%
  na.omit
CG_4=CG_4 %>%
  na_if("") %>%
  na.omit
#Drop contig
CG_1=subset(CG_1,select=-c(contig))
CG_2=subset(CG_2,select=-c(contig))
CG_3=subset(CG_3,select=-c(contig))
CG_4=subset(CG_4,select=-c(contig))

#Make new df
rpk1=CG_1
rpk2=CG_2
rpk3=CG_3
rpk4=CG_4
#Put length in KB
rpk1$Length=rpk1$Length/1000
rpk2$Length=rpk2$Length/1000
rpk3$Length=rpk3$Length/1000
rpk4$Length=rpk4$Length/1000
#Create RPK for Mat 1
rpk1$R1_1=rpk1$R1_1/rpk1$Length
rpk1$R1_2=rpk1$R1_2/rpk1$Length
rpk1$R1_3=rpk1$R1_3/rpk1$Length
rpk1$R1_4=rpk1$R1_4/rpk1$Length
rpk1$R1_5=rpk1$R1_5/rpk1$Length
#Create RPK for Mat 2
rpk2$R2_1=rpk2$R2_1/rpk2$Length
rpk2$R2_2=rpk2$R2_2/rpk2$Length
rpk2$R2_3=rpk2$R2_3/rpk2$Length
rpk2$R2_4=rpk2$R2_4/rpk2$Length
rpk2$R2_5=rpk2$R2_5/rpk2$Length
#Create RPK for Mat 3
rpk3$R3_1=rpk3$R3_1/rpk3$Length
rpk3$R3_2=rpk3$R3_2/rpk3$Length
rpk3$R3_3=rpk3$R3_3/rpk3$Length
rpk3$R3_4=rpk3$R3_4/rpk3$Length
rpk3$R3_5=rpk3$R3_5/rpk3$Length
#Create RPK for Mat 4
rpk4$R4_1=rpk4$R4_1/rpk4$Length
rpk4$R4_2=rpk4$R4_2/rpk4$Length
rpk4$R4_3=rpk4$R4_3/rpk4$Length
rpk4$R4_4=rpk4$R4_4/rpk4$Length
rpk4$R4_5=rpk4$R4_5/rpk4$Length
#create scaling factors
#Mat 1
sum11=(sum(rpk1$R1_1)/1000000)
sum12=(sum(rpk1$R1_2)/1000000)
sum13=(sum(rpk1$R1_3)/1000000)
sum14=(sum(rpk1$R1_4)/1000000)
sum15=(sum(rpk1$R1_5)/1000000)
#Mat 2
sum21=(sum(rpk2$R2_1)/1000000)
sum22=(sum(rpk2$R2_2)/1000000)
sum23=(sum(rpk2$R2_3)/1000000)
sum24=(sum(rpk2$R2_4)/1000000)
sum25=(sum(rpk2$R2_5)/1000000)
#Mat 3
sum31=(sum(rpk3$R3_1)/1000000)
sum32=(sum(rpk3$R3_2)/1000000)
sum33=(sum(rpk3$R3_3)/1000000)
sum34=(sum(rpk3$R3_4)/1000000)
sum35=(sum(rpk3$R3_5)/1000000)
#Mat 4
sum41=(sum(rpk4$R4_1)/1000000)
sum42=(sum(rpk4$R4_2)/1000000)
sum43=(sum(rpk4$R4_3)/1000000)
sum44=(sum(rpk4$R4_4)/1000000)
sum45=(sum(rpk4$R4_5)/1000000)
#Create TPM for Mat 1
rpk1$R1_1=rpk1$R1_1/sum11
rpk1$R1_2=rpk1$R1_2/sum12
rpk1$R1_3=rpk1$R1_3/sum13
rpk1$R1_4=rpk1$R1_4/sum14
rpk1$R1_5=rpk1$R1_5/sum15
#Create TPM for Mat 2
rpk2$R2_1=rpk2$R2_1/sum21
rpk2$R2_2=rpk2$R2_2/sum22
rpk2$R2_3=rpk2$R2_3/sum23
rpk2$R2_4=rpk2$R2_4/sum24
rpk2$R2_5=rpk2$R2_5/sum25
#Create TPM for Mat 3
rpk3$R3_1=rpk3$R3_1/sum31
rpk3$R3_2=rpk3$R3_2/sum32
rpk3$R3_3=rpk3$R3_3/sum33
rpk3$R3_4=rpk3$R3_4/sum34
rpk3$R3_5=rpk3$R3_5/sum35
#Create TPM for Mat 4
rpk4$R4_1=rpk4$R4_1/sum41
rpk4$R4_2=rpk4$R4_2/sum42
rpk4$R4_3=rpk4$R4_3/sum43
rpk4$R4_4=rpk4$R4_4/sum44
rpk4$R4_5=rpk4$R4_5/sum45

tpm1=rpk1 %>%
  group_by(kegg_id,Class) %>%
  summarise(across(1:5, sum))
tpm2=rpk2 %>%
  group_by(kegg_id,Class) %>%
  summarise(across(1:5, sum))
tpm3=rpk3 %>%
  group_by(kegg_id,Class) %>%
  summarise(across(1:5, sum))
tpm4=rpk4 %>%
  group_by(kegg_id,Class) %>%
  summarise(across(1:5, sum))

#Split into 2 dataframes each with 1) genes of interest for normalization and 2) everything else
norml1=tpm1[which(tpm1$kegg_id=="K03553" | tpm1$kegg_id=="K03043" | tpm1$kegg_id=="K03046" | tpm1$kegg_id=="K03070"),]
norml1=as.data.frame(norml1)
norml2=tpm2[which(tpm2$kegg_id=="K03553" | tpm2$kegg_id=="K03043" | tpm2$kegg_id=="K03046" | tpm2$kegg_id=="K03070"),]
norml2=as.data.frame(norml2)
norml3=tpm3[which(tpm3$kegg_id=="K03553" | tpm3$kegg_id=="K03043" | tpm3$kegg_id=="K03046" | tpm3$kegg_id=="K03070"),]
norml3=as.data.frame(norml3)
norml4=tpm4[which(tpm4$kegg_id=="K03553" | tpm4$kegg_id=="K03043" | tpm4$kegg_id=="K03046" | tpm4$kegg_id=="K03070"),]
norml4=as.data.frame(norml4)

othr1=tpm1[which(tpm1$kegg_id!="K03553" | tpm1$kegg_id!="K03043" | tpm1$kegg_id!="K03046" | tpm1$kegg_id!="K03070"),]
othr1=as.data.frame(othr1)
othr2=tpm2[which(tpm2$kegg_id!="K03553" | tpm2$kegg_id!="K03043" | tpm2$kegg_id!="K03046" | tpm2$kegg_id!="K03070"),]
othr2=as.data.frame(othr2)
othr3=tpm3[which(tpm3$kegg_id!="K03553" | tpm3$kegg_id!="K03043" | tpm3$kegg_id!="K03046" | tpm3$kegg_id!="K03070"),]
othr3=as.data.frame(othr3)
othr4=tpm4[which(tpm4$kegg_id!="K03553" | tpm4$kegg_id!="K03043" | tpm4$kegg_id!="K03046" | tpm4$kegg_id!="K03070"),]
othr4=as.data.frame(othr4)

#Find mean expression for each Class
mn1=norml1%>%
  group_by(Class) %>%
  summarise_at(vars("R1_1","R1_2","R1_3","R1_4","R1_5"),mean)
mn2=norml2%>%
  group_by(Class) %>%
  summarise_at(vars("R2_1","R2_2","R2_3","R2_4","R2_5"),mean)

#Merge to get metabolic subset
c1=merge(energy1,othr1, by="kegg_id")
c2=merge(energy2,othr2,by="kegg_id")
c3=merge(energy3,othr3, by="kegg_id")
c4=merge(energy4,othr4, by="kegg_id")

#Merge
CG_merge1=merge(y=c1,x=c2,by=c("kegg_id","Class"),all=T)
CG_merge2=merge(c3,c4,by=c("kegg_id","Class"),all=T)
CG_merge=merge(y=CG_merge1,x=CG_merge2,by=c("kegg_id","Class"),all=T)

#Replace NA with 0
CG_merge[is.na(CG_merge)] <- 0

#Convert to long
MT_long <- gather(MT_merge, Sample, TPM, R3_1:R1_5)

#Add column for time and MatID
conflict_prefer("mutate", "dplyr")
conflict_prefer("startsWith", "gdata")

MT_long=MT_long %>%
  mutate(Time=case_when(
    endsWith(Sample,"1") ~ "1",
    endsWith(Sample,"2") ~ "2",
    endsWith(Sample,"3") ~ "3",
    endsWith(Sample,"4") ~ "4",
    endsWith(Sample,"5") ~ "5",
  ))
MT_long=MT_long %>%
  mutate(Mat=case_when(
    str_detect(Sample,"R1") ~ "1",
    str_detect(Sample,"R2") ~ "2",
    str_detect(Sample,"R3") ~ "3",
    str_detect(Sample,"R4") ~ "4",
    str_detect(Sample,"R5") ~ "5",
  ))








############
#TRY RNA WITHOUT FILTERING BY KNOWN KEGG
###########
#####DIEL RNA DATA ALIGNED TO MAGS
#We will have to do some merging to get the feature count matrix
#Read in all the counts
rna_diel1=read.xls('/Users/cissell/Desktop/D1_RNA_counts.xlsx')
rna_diel2=read.xls('/Users/cissell/Desktop/D2_RNA_counts.xlsx')
rna_diel3=read.xls('/Users/cissell/Desktop/D3_RNA_counts.xlsx')
rna_diel4=read.xls('/Users/cissell/Desktop/D4_RNA_counts.xlsx')
#Read in the gene id
rna_kegg1=read.xls('/Users/cissell/Desktop/D1_annotations.xlsx')
rna_kegg2=read.xls('/Users/cissell/Desktop/D2_annotations.xlsx')
rna_kegg3=read.xls('/Users/cissell/Desktop/D3_annotations.xlsx')
rna_kegg4=read.xls('/Users/cissell/Desktop/D4_annotations.xlsx')
#Merge counts with gene ids
merged1=merge(rna_diel1,rna_kegg1, by="contig",all.x=T)
merged2=merge(rna_diel2,rna_kegg2, by="contig",all.x=T)
merged3=merge(rna_diel3,rna_kegg3, by="contig",all.x=T)
merged4=merge(rna_diel4,rna_kegg4, by="contig",all.x=T)
#Remove rows without KEGG ID
#Drop contig
merged1=subset(merged1,select=-c(contig))
merged2=subset(merged2,select=-c(contig))
merged3=subset(merged3,select=-c(contig))
merged4=subset(merged4,select=-c(contig))
#Drop length
merged1=subset(merged1,select=-c(Length))
merged2=subset(merged2,select=-c(Length))
merged3=subset(merged3,select=-c(Length))
merged4=subset(merged4,select=-c(Length))
#Replace Blank KEGG
merged1$kegg_id = sub("^$","Unknown",merged1$kegg_id)
merged2$kegg_id = sub("^$","Unknown",merged2$kegg_id)
merged3$kegg_id = sub("^$","Unknown",merged3$kegg_id)
merged4$kegg_id = sub("^$","Unknown",merged4$kegg_id)
#Summarize within KEGG ID
summerged1=merged1 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
summerged2=merged2 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
summerged3=merged3 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
summerged4=merged4 %>%
  group_by(kegg_id) %>%
  summarise(across(1:5, sum))
#Merge across gene ids
rna_merge1=merge(y=summerged1,x=summerged2,by="kegg_id",all.x=T)
rna_merge2=merge(summerged3,summerged4,by="kegg_id",all.x=T)
rna_merge=merge(y=rna_merge1,x=rna_merge2,by="kegg_id",all.x=T)

rna_merge[is.na(rna_merge)]=0
#Make OTU table for phyloseq
#Rownames have to be ID
rownames(rna_merge)<- rna_merge$kegg_id

#Check counts
count_rna <- colSums(rna_merge[, c(2:ncol(rna_merge))])
count_rna

#Create taxonomy table
rna_tax<- rna_merge %>%
  select(kegg_id)
str(rna_tax)

#Make factors not strings
rna_tax<- rna_tax %>%
  mutate_if(is.character, as.factor)
str(rna_tax)

#Rename first column
colnames(rna_tax)[1]<- "kegg_id"
str(rna_tax)
#OTU IDs have to be rownames, just like OTU table
rownames(rna_tax)<- rna_tax$kegg_id

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
rna_merge<- subset(rna_merge, select=-kegg_id)


#Create metadata table
rna_meta=read.xls('/Users/cissell/Desktop/rna_diel_meta.xlsx')
str(rna_meta)
rownames(rna_meta)<- rna_meta$SampleID
rna_meta<- rna_meta %>%
  select(-SampleID)
#Make proper levels
rna_meta$Mat<- factor(rna_meta$Mat, levels = c("R1", "R2", "R3","R4"))
rna_meta$Age<- factor(rna_meta$Age, levels = c("1", "2","3","4", "5"))

#Create matrices before creating phyloseq object
rna_otu_mat<- as.matrix(rna_merge)
rna_tax_mat<- tax_table(as.matrix(rna_tax))

#Transform to phyloseq objects
rna_phylo_OTU<- otu_table(rna_otu_mat, taxa_are_rows = TRUE)
rna_phylo_TAX<- tax_table(rna_tax_mat)
rna_phylo_samples<- sample_data(rna_meta)

#And unite them into one object
rna_phylo_object<- phyloseq(rna_phylo_OTU, rna_phylo_TAX, rna_phylo_samples)
#Check and see if it makes sense
sample_sums(rna_phylo_object)
sample_names(rna_phylo_object)
rank_names(rna_phylo_object)
sample_variables(rna_phylo_object)
otu_table(rna_phylo_object)[1:3, 1:2]
taxa_names(rna_phylo_object)[1:5]
tax_table(rna_phylo_object)
taxa_names(rna_phylo_object)

library(microbiome)
rna_clr=microbiome::transform(rna_phylo_object, transform='clr',target="sample")
#NMDS on raw dataframe
#MAke dataframe
rna_nmds_df=psmelt(rna_clr)
rna_nmds_df=rna_nmds_df%>%
  select(-OTU)%>%
  spread(rna_nmds_df,key=kegg_id,value=Abundance)

dat=rna_nmds_df[,4:6001]
dat=as.data.frame(lapply(dat,unlist))
str(dat)
grp=as.data.frame(rna_nmds_df[,2:3],header=FALSE)
rna_nmds=metaMDS(dat,distance="euclidean")
rna_nmds
#Stressplot
stressplot(rna_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)
dev.print(png,fil="stress_rna_full.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
#Scree plot of stress
library(goeveg)
dimcheckMDS(
  dat,
  distance = "euclidean",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

kegg.fit=envfit(rna_nmds,dat,permutations=999)
plot(rna_nmds,type="t")
data.score=as.data.frame(scores(rna_nmds))
data.score$site=rownames(data.score)
data.score$grp=grp

species.scores=envfit(rna_nmds,dat)
species.scores2=data.frame((species.scores$vectors)$arrows,(species.scores$vectors)$r,(species.scores$vectors)$pvals)
species.scores2=species.scores2[,-c(3)]
names(species.scores2)[names(species.scores2)=="X.species.scores.vectors..pvals"]="pval"
head(species.scores2)
sig.sp.sc=subset(species.scores2,pval<=0.05)
sig.sp.sc=rownames_to_column(sig.sp.sc, "KEGG")


#Make hull vectors and hull data frame for plotting polygons
age.1 = data.score[(data.score$grp)$Age == "1",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.2 = data.score[(data.score$grp)$Age == "2",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.3=data.score[(data.score$grp)$Age == "3",][chull(data.score[(data.score$grp)$Age== "3",c("NMDS1","NMDS2")]),]
age.4=data.score[(data.score$grp)$Age == "4",][chull(data.score[(data.score$grp)$Age== "4",c("NMDS1","NMDS2")]),]
age.5=data.score[(data.score$grp)$Age == "5",][chull(data.score[(data.score$grp)$Age=="5",c("NMDS1","NMDS2")]),]
hull.data=rbind(age.1,age.2,age.3,age.4,age.5)
str(hull.data)


mat.1 = data.score[(data.score$grp)$Mat == "R1",][chull(data.score[(data.score$grp)$Mat=="R1",c("NMDS1","NMDS2")]),]
mat.2 = data.score[(data.score$grp)$Mat == "R2",][chull(data.score[(data.score$grp)$Mat=="R2",c("NMDS1","NMDS2")]),]
mat.3=data.score[(data.score$grp)$Mat == "R3",][chull(data.score[(data.score$grp)$Mat== "R3",c("NMDS1","NMDS2")]),]
mat.4=data.score[(data.score$grp)$Mat == "R4",][chull(data.score[(data.score$grp)$Mat== "R4",c("NMDS1","NMDS2")]),]
mat.data=rbind(mat.1,mat.2,mat.3,mat.4)
str(mat.data)
#Plot
nmds_rna=ggplot() +
  #geom_polygon(data=hull.data,
  #             aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
  #             alpha=0.4)+
  geom_polygon(data=mat.data,
               aes(x=NMDS1,y=NMDS2, group=grp$Mat),
               #size=1,
               #linetype=1,
               fill="#333333",
               #colour="#333333",
               alpha=0.6)+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2"),name="Time") +
  #  scale_fill_manual(values=c("#90a295","#fde992","#455765","#add8e6","#cd9fb2","#fde992")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=49,y=-37, label="Stress = 0.164", color="black",size=5)+
  #annotate(geom="text",x=0,y=-1, label="3", color="#333333",size=5)+
  #annotate(geom="text",x=0.78,y=-0.42, label="1", color="#333333",size=5)+
  #annotate(geom="text",x=-0.8,y=0.2, label="4", color="#333333",size=5)+
  #annotate(geom="text",x=0.3,y=1.25, label="2", color="#333333",size=5)+
  scale_x_continuous(name="NMDS1")+
  scale_y_continuous(name="NMDS2")
nmds_rna

dev.print(png,fil="nmds_rna_full.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

#Create distance matrix from clr
clr_dist_matrix <- phyloseq::distance(rna_clr, method = "euclidean")
distance(phylo_clr, method = "euclidean")

#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Age)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.8047

#Mat
dispr2 <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(rna_clr)$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.9847
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(clr_dist_matrix ~sample_data(rna_clr)$Age+sample_data(rna_clr)$Mat,permutations=20000, method="euclidean", by="margin")
permanova

#Pairwise PERMANOVA
#Time
pair_perm=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Age, sim.method="euclidean")
pair_perm #none
#Mat
pair_perm2=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Mat, sim.method="euclidean")
pair_perm2 #2v4 and 3v4





########################
#################
#######
#PHAGES RNA
#######
#################
########################
#We will have to do some merging to get the feature count matrix
#Read in all the counts
phage_diel1=read.xls('/Users/cissell/Desktop/D1_RNA_phages_counts.xlsx')
phage_diel2=read.xls('/Users/cissell/Desktop/D2_RNA_phages_counts.xlsx')
phage_diel3=read.xls('/Users/cissell/Desktop/D3_RNA_phages_counts.xlsx')
phage_diel4=read.xls('/Users/cissell/Desktop/D4_RNA_phages_counts.xlsx')
#Read in the gene id
phage_kegg1=read.xls('/Users/cissell/Desktop/D1_phage_annotations.xlsx')
phage_kegg2=read.xls('/Users/cissell/Desktop/D2_phage_annotations.xlsx')
phage_kegg3=read.xls('/Users/cissell/Desktop/D3_phage_annotations.xlsx')
phage_kegg4=read.xls('/Users/cissell/Desktop/D4_phage_annotations.xlsx')
#Merge counts with gene ids
merged1=merge(phage_diel1,phage_kegg1, by="contig",all.x=T)
merged2=merge(phage_diel2,phage_kegg2, by="contig",all.x=T)
merged3=merge(phage_diel3,phage_kegg3, by="contig",all.x=T)
merged4=merge(phage_diel4,phage_kegg4, by="contig",all.x=T)
#Remove rows without KEGG ID
smerged1=merged1 %>%
  na_if("") %>%
  na.omit
smerged2=merged2 %>%
  na_if("") %>%
  na.omit
smerged3=merged3 %>%
  na_if("") %>%
  na.omit
smerged4=merged4 %>%
  na_if("") %>%
  na.omit
#Drop contig
smerged1=subset(smerged1,select=-c(contig,Length,viral_hit))
smerged2=subset(smerged2,select=-c(contig,Length,viral_hit))
smerged3=subset(smerged3,select=-c(contig,Length,viral_hit))
smerged4=subset(smerged4,select=-c(contig,Length,viral_hit))
#Summarize within KEGG ID
summerged1=smerged1 %>%
  group_by(viral_id) %>%
  summarise(across(1:5, sum))
summerged2=smerged2 %>%
  group_by(viral_id) %>%
  summarise(across(1:5, sum))
summerged3=smerged3 %>%
  group_by(viral_id) %>%
  summarise(across(1:5, sum))
summerged4=smerged4 %>%
  group_by(viral_id) %>%
  summarise(across(1:5, sum))
#Merge across gene ids
rna_merge1=merge(x=summerged1,y=summerged2,by="viral_id")
rna_merge2=merge(y=summerged3,x=summerged4,by="viral_id")
rna_merge=merge(x=rna_merge1,y=rna_merge2,by="viral_id")

#Replace NA with 0
rna_merge[is.na(rna_merge)]=0

#Make OTU table for phyloseq
#Rownames have to be ID
rownames(rna_merge)<- rna_merge$viral_id

#Check counts
count_rna <- colSums(rna_merge[, c(2:ncol(rna_merge))])
count_rna

#Create taxonomy table
rna_tax<- rna_merge %>%
  select(viral_id)
str(rna_tax)

#Make factors not strings
rna_tax<- rna_tax %>%
  mutate_if(is.character, as.factor)
str(rna_tax)

#Rename first column
colnames(rna_tax)[1]<- "viral_id"
str(rna_tax)
#OTU IDs have to be rownames, just like OTU table
rownames(rna_tax)<- rna_tax$viral_id

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
rna_merge<- subset(rna_merge, select=-viral_id)


#Create metadata table
rna_meta=read.xls('/Users/cissell/Desktop/rna_diel_meta.xlsx')
str(rna_meta)
rownames(rna_meta)<- rna_meta$SampleID
rna_meta<- rna_meta %>%
  select(-SampleID)
#Make proper levels
rna_meta$Mat<- factor(rna_meta$Mat, levels = c("R1", "R2", "R3","R4"))
rna_meta$Age<- factor(rna_meta$Age, levels = c("1", "2","3","4", "5"))

#Create matrices before creating phyloseq object
rna_otu_mat<- as.matrix(rna_merge)
rna_tax_mat<- tax_table(as.matrix(rna_tax))

#Transform to phyloseq objects
rna_phylo_OTU<- otu_table(rna_otu_mat, taxa_are_rows = TRUE)
rna_phylo_TAX<- tax_table(rna_tax_mat)
rna_phylo_samples<- sample_data(rna_meta)

#And unite them into one object
rna_phylo_object<- phyloseq(rna_phylo_OTU, rna_phylo_TAX, rna_phylo_samples)
#Check and see if it makes sense
sample_sums(rna_phylo_object)
sample_names(rna_phylo_object)
rank_names(rna_phylo_object)
sample_variables(rna_phylo_object)
otu_table(rna_phylo_object)[1:3, 1:2]
taxa_names(rna_phylo_object)[1:5]
tax_table(rna_phylo_object)
taxa_names(rna_phylo_object)

library(microbiome)
phage_rna_clr=microbiome::transform(rna_phylo_object, transform='clr',target="sample")
#NMDS on raw dataframe
#MAke dataframe
rna_nmds_df=psmelt(phage_rna_clr)
rna_nmds_df=rna_nmds_df%>%
  select(-OTU)%>%
  spread(rna_nmds_df,key=viral_id,value=Abundance)

dat=rna_nmds_df[,4:466]
dat=as.data.frame(lapply(dat,unlist))
str(dat)
grp=as.data.frame(rna_nmds_df[,2:3],header=FALSE)
rna_nmds=metaMDS(dat,distance="euclidean")
rna_nmds
#Stressplot
stressplot(rna_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)
dev.print(png,fil="stress_rna_full.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
#Scree plot of stress
library(goeveg)
dimcheckMDS(
  dat,
  distance = "euclidean",
  k = 6,
  trymax = 20,
  autotransform = TRUE
)

kegg.fit=envfit(rna_nmds,dat,permutations=999)
plot(rna_nmds,type="t")
data.score=as.data.frame(scores(rna_nmds))
data.score$site=rownames(data.score)
data.score$grp=grp

species.scores=envfit(rna_nmds,dat)
species.scores2=data.frame((species.scores$vectors)$arrows,(species.scores$vectors)$r,(species.scores$vectors)$pvals)
species.scores2=species.scores2[,-c(3)]
names(species.scores2)[names(species.scores2)=="X.species.scores.vectors..pvals"]="pval"
head(species.scores2)
sig.sp.sc=subset(species.scores2,pval<=0.05)
sig.sp.sc=rownames_to_column(sig.sp.sc, "Vir_ID")


#Make hull vectors and hull data frame for plotting polygons
age.1 = data.score[(data.score$grp)$Age == "1",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.2 = data.score[(data.score$grp)$Age == "2",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.3=data.score[(data.score$grp)$Age == "3",][chull(data.score[(data.score$grp)$Age== "3",c("NMDS1","NMDS2")]),]
age.4=data.score[(data.score$grp)$Age == "4",][chull(data.score[(data.score$grp)$Age== "4",c("NMDS1","NMDS2")]),]
age.5=data.score[(data.score$grp)$Age == "5",][chull(data.score[(data.score$grp)$Age=="5",c("NMDS1","NMDS2")]),]
hull.data=rbind(age.1,age.2,age.3,age.4,age.5)
str(hull.data)


mat.1 = data.score[(data.score$grp)$Mat == "R1",][chull(data.score[(data.score$grp)$Mat=="R1",c("NMDS1","NMDS2")]),]
mat.2 = data.score[(data.score$grp)$Mat == "R2",][chull(data.score[(data.score$grp)$Mat=="R2",c("NMDS1","NMDS2")]),]
mat.3=data.score[(data.score$grp)$Mat == "R3",][chull(data.score[(data.score$grp)$Mat== "R3",c("NMDS1","NMDS2")]),]
mat.4=data.score[(data.score$grp)$Mat == "R4",][chull(data.score[(data.score$grp)$Mat== "R4",c("NMDS1","NMDS2")]),]
mat.data=rbind(mat.1,mat.2,mat.3,mat.4)
str(mat.data)
#Plot
nmds_rna=ggplot() +
  #geom_polygon(data=hull.data,
  #             aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
  #             alpha=0.4)+
  geom_polygon(data=mat.data,
               aes(x=NMDS1,y=NMDS2, group=grp$Mat),
               #size=1,
               #linetype=1,
               fill="#333333",
               #colour="#333333",
               alpha=0.6)+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2"),name="Time") +
  #  scale_fill_manual(values=c("#90a295","#fde992","#455765","#add8e6","#cd9fb2","#fde992")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=49,y=-37, label="Stress = 0.164", color="black",size=5)+
  #annotate(geom="text",x=0,y=-1, label="3", color="#333333",size=5)+
  #annotate(geom="text",x=0.78,y=-0.42, label="1", color="#333333",size=5)+
  #annotate(geom="text",x=-0.8,y=0.2, label="4", color="#333333",size=5)+
  #annotate(geom="text",x=0.3,y=1.25, label="2", color="#333333",size=5)+
  scale_x_continuous(name="NMDS1")+
  scale_y_continuous(name="NMDS2")
nmds_rna

dev.print(png,fil="nmds_rna_full.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

#Create distance matrix from clr
phage_dist_matrix <- phyloseq::distance(phage_rna_clr, method = "euclidean")
distance(phylo_clr, method = "euclidean")

#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(phage_dist_matrix, phyloseq::sample_data(rna_clr)$Age)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.8047

#Mat
dispr2 <- vegan::betadisper(phage_dist_matrix, phyloseq::sample_data(rna_clr)$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.9847
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(phage_dist_matrix ~sample_data(rna_clr)$Age+sample_data(rna_clr)$Mat,permutations=999, method="euclidean", by="margin")
permanova

#Pairwise PERMANOVA
#Time
pair_perm=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Age, sim.method="euclidean")
pair_perm #none
#Mat
pair_perm2=pairwiseAdonis::pairwise.adonis(clr_dist_matrix, factors=sample_data(rna_clr)$Mat, sim.method="euclidean")
pair_perm2 #2v4 and 3v4


library(MicEco)
#Time
venn_age=ps_venn(rna_phylo_object,group='Age',type='counts',relative=FALSE,fraction=0,fill=c("#90a295","#fde992","#455765","#add8e6","#cd9fb2","#fde992"))
#Individual
venn_mat=ps_venn(rna_phylo_object,group='Mat',type='percent',relative=FALSE,fraction=0,fill=c("#90a295","#455765","#cd9fb2","#fde992"))
#Combine
venns=ggarrange(venn_age,venn_mat,
                ncol=2,nrow=1,align="hv",
                labels=c("a","b"))
venns
annotate_figure(venns,
                left=text_grob("Phylum", color="Black", rot=90,size=15,vjust=(2.3),hjust=0.1))
dev.print(png,fil="venns.png",type="quartz",antialias="default",width=6.5,
          height=4.5,units="in",res=1300)




#########################
#####RNA PHAGES LYSOGENY HALLMARK MARKERS
########################
phage_diel1=read.xls('/Users/cissell/Desktop/D1_RNA_phages_counts.xlsx')
phage_diel2=read.xls('/Users/cissell/Desktop/D2_RNA_phages_counts.xlsx')
phage_diel3=read.xls('/Users/cissell/Desktop/D3_RNA_phages_counts.xlsx')
phage_diel4=read.xls('/Users/cissell/Desktop/D4_RNA_phages_counts.xlsx')
#Read in the gene id
phage_kegg1=read.xls('/Users/cissell/Desktop/D1_phage_annotations.xlsx')
phage_kegg2=read.xls('/Users/cissell/Desktop/D2_phage_annotations.xlsx')
phage_kegg3=read.xls('/Users/cissell/Desktop/D3_phage_annotations.xlsx')
phage_kegg4=read.xls('/Users/cissell/Desktop/D4_phage_annotations.xlsx')
#Merge counts with gene ids
merged1=merge(phage_diel1,phage_kegg1, by="contig",all.x=T)
merged2=merge(phage_diel2,phage_kegg2, by="contig",all.x=T)
merged3=merge(phage_diel3,phage_kegg3, by="contig",all.x=T)
merged4=merge(phage_diel4,phage_kegg4, by="contig",all.x=T)
#Remove rows without KEGG ID
smerged1=merged1 %>%
  na_if("") %>%
  na.omit
smerged2=merged2 %>%
  na_if("") %>%
  na.omit
smerged3=merged3 %>%
  na_if("") %>%
  na.omit
smerged4=merged4 %>%
  na_if("") %>%
  na.omit
#Separate phage to combine with host-linkages
sep1=smerged1%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11"))
sep1=merge(sep1,D1_ph_count_abund_merge,by="Phage")
sep1=sep1%>%
  mutate_at("Phylum",str_replace,"p__","")

sep2=smerged2%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11"))
sep2=merge(sep2,D2_ph_indiv_abund_merge,by="Phage")
sep2=sep2%>%
  mutate_at("Phylum",str_replace,"p__","")

sep3=smerged3%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11"))
sep3=merge(sep3,D3_ph_indiv_abund_merge,by="Phage")
sep3=sep3%>%
  mutate_at("Phylum",str_replace,"p__","")

sep4=smerged4%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11"))
sep4=merge(sep4,D4_ph_indiv_abund_merge,by="Phage")
sep4=sep4%>%
  mutate_at("Phylum",str_replace,"p__","")
#Now drop unecessary columns added when merging
#sep1=sep1%>%
#  select(-c("H1_1","H1_3","H1_5","P1_1","P1_3","P1_5","Host"))
#sep2=sep2%>%
#  select(-c("H2_1","H2_3","H2_5","P2_1","P2_3","P2_5","Host"))
#sep3=sep3%>%
#  select(-c("H3_1","H3_3","H3_5","P3_1","P3_3","P3_5","Host"))
#sep4=sep4%>%
#  select(-c("H4_1","H4_3","H4_5","P4_1","P4_3","P4_5","Host"))

#Select for integrase, recombinase  and ParA/ParB genes (lysogeny)
sep1_lysogen=sep1[str_detect(sep1$viral_hit,
                                  "integrase | ParA family | ParB | recombinase"),]
sep2_lysogen=sep2[str_detect(sep2$viral_hit,
                             "integrase | ParA family | ParB | recombinase"),]
sep3_lysogen=sep3[str_detect(sep3$viral_hit,
                             "integrase | ParA family | ParB | recombinase"),]
sep4_lysogen=sep4[str_detect(sep4$viral_hit,
                             "integrase | ParA family | ParB | recombinase"),]

#For lysis, lets look at Terminase which is involved in packing procapsids
sep1_lysis=sep1[str_detect(sep1$viral_hit,
                            "terminase"),]
sep2_lysis=sep2[str_detect(sep2$viral_hit,
                           "terminase"),]
sep3_lysis=sep3[str_detect(sep3$viral_hit,
                           "terminase"),]
sep4_lysis=sep4[str_detect(sep4$viral_hit,
                           "terminase"),]

#Summarize by phyla to get feature called "Lysis" and "Lysogeny" which represents summed expression of all lysogeny and lysis genes per phyla for proportionality analysis
#D1_lysogeny
sep1_lysogen=sep1_lysogen%>%
  select(R1_1,R1_3,R1_5,P1_1,P1_3,P1_5,H1_1,H1_3,H1_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep1_lysogen=pivot_longer(data=sep1_lysogen,
               cols=1:9,
               names_to=c("Type","Time"),
               names_sep="_")

sep1_lysogen=pivot_wider(data=sep1_lysogen,
                         names_from=Type,
                         values_from=value)

sep1_lysogen=sep1_lysogen%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D1_1")%>%
  mutate_at("Sample",str_replace,"3","D1_3")%>%
  mutate_at("Sample",str_replace,"5","D1_5")%>%
  rename(Lysogeny=R1,
         Phage=P1,
         Host=H1)
#D1_lysis
sep1_lysis=sep1_lysis%>%
  select(R1_1,R1_3,R1_5,P1_1,P1_3,P1_5,H1_1,H1_3,H1_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep1_lysis=pivot_longer(data=sep1_lysis,
                          cols=1:9,
                          names_to=c("Type","Time"),
                          names_sep="_")

sep1_lysis=pivot_wider(data=sep1_lysis,
                         names_from=Type,
                         values_from=value)

sep1_lysis=sep1_lysis%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D1_1")%>%
  mutate_at("Sample",str_replace,"3","D1_3")%>%
  mutate_at("Sample",str_replace,"5","D1_5")%>%
  rename(Lysis=R1,
         Phage=P1,
         Host=H1)%>%
  select(Lysis)
#Combine lysis into lysogeny for mat 1
sep1_ll=cbind(sep1_lysogen,sep1_lysis)

#Now mat 2
#D2_lysogeny
sep2_lysogen=sep2_lysogen%>%
  select(R2_1,R2_3,R2_5,P2_1,P2_3,P2_5,H2_1,H2_3,H2_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep2_lysogen=pivot_longer(data=sep2_lysogen,
                          cols=1:9,
                          names_to=c("Type","Time"),
                          names_sep="_")

sep2_lysogen=pivot_wider(data=sep2_lysogen,
                         names_from=Type,
                         values_from=value)

sep2_lysogen=sep2_lysogen%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D2_1")%>%
  mutate_at("Sample",str_replace,"3","D2_3")%>%
  mutate_at("Sample",str_replace,"5","D2_5")%>%
  rename(Lysogeny=R2,
         Phage=P2,
         Host=H2)
#D2_lysis
sep2_lysis=sep2_lysis%>%
  select(R2_1,R2_3,R2_5,P2_1,P2_3,P2_5,H2_1,H2_3,H2_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep2_lysis=pivot_longer(data=sep2_lysis,
                        cols=1:9,
                        names_to=c("Type","Time"),
                        names_sep="_")

sep2_lysis=pivot_wider(data=sep2_lysis,
                       names_from=Type,
                       values_from=value)

sep2_lysis=sep2_lysis%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D2_1")%>%
  mutate_at("Sample",str_replace,"3","D2_3")%>%
  mutate_at("Sample",str_replace,"5","D2_5")%>%
  rename(Lysis=R2,
         Phage=P2,
         Host=H2)%>%
  select(Lysis)
#Combine lysis into lysogeny for mat 2
sep2_ll=cbind(sep2_lysogen,sep2_lysis)


#Now mat 3
#D3_lysogeny
sep3_lysogen=sep3_lysogen%>%
  select(R3_1,R3_3,R3_5,P3_1,P3_3,P3_5,H3_1,H3_3,H3_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep3_lysogen=pivot_longer(data=sep3_lysogen,
                          cols=1:9,
                          names_to=c("Type","Time"),
                          names_sep="_")

sep3_lysogen=pivot_wider(data=sep3_lysogen,
                         names_from=Type,
                         values_from=value)

sep3_lysogen=sep3_lysogen%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D3_1")%>%
  mutate_at("Sample",str_replace,"3","D3_3")%>%
  mutate_at("Sample",str_replace,"5","D3_5")%>%
  mutate_at("Sample",str_replace,"DD3_3_1","D3_1")%>%
  rename(Lysogeny=R3,
         Phage=P3,
         Host=H3)
#D3_lysis
sep3_lysis=sep3_lysis%>%
  select(R3_1,R3_3,R3_5,P3_1,P3_3,P3_5,H3_1,H3_3,H3_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep3_lysis=pivot_longer(data=sep3_lysis,
                        cols=1:9,
                        names_to=c("Type","Time"),
                        names_sep="_")

sep3_lysis=pivot_wider(data=sep3_lysis,
                       names_from=Type,
                       values_from=value)

sep3_lysis=sep3_lysis%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D3_1")%>%
  mutate_at("Sample",str_replace,"3","D3_3")%>%
  mutate_at("Sample",str_replace,"5","D3_5")%>%
  rename(Lysis=R3,
         Phage=P3,
         Host=H3)%>%
  select(Lysis)
#Combine lysis into lysogeny for mat 3
sep3_ll=cbind(sep3_lysogen,sep3_lysis)


#Now mat 4
#D4_lysogeny
sep4_lysogen=sep4_lysogen%>%
  select(R4_1,R4_3,R4_5,P4_1,P4_3,P4_5,H4_1,H4_3,H4_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep4_lysogen=pivot_longer(data=sep4_lysogen,
                          cols=1:9,
                          names_to=c("Type","Time"),
                          names_sep="_")

sep4_lysogen=pivot_wider(data=sep4_lysogen,
                         names_from=Type,
                         values_from=value)

sep4_lysogen=sep4_lysogen%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D4_1")%>%
  mutate_at("Sample",str_replace,"3","D4_3")%>%
  mutate_at("Sample",str_replace,"5","D4_5")%>%
  rename(Lysogeny=R4,
         Phage=P4,
         Host=H4)
#D4_lysis
sep4_lysis=sep4_lysis%>%
  select(R4_1,R4_3,R4_5,P4_1,P4_3,P4_5,H4_1,H4_3,H4_5,Phylum)%>%
  group_by(Phylum)%>%
  summarise_all(sum)%>%
  ungroup()%>%
  select(-Phylum)%>%
  summarise_all(sum)

sep4_lysis=pivot_longer(data=sep4_lysis,
                        cols=1:9,
                        names_to=c("Type","Time"),
                        names_sep="_")

sep4_lysis=pivot_wider(data=sep4_lysis,
                       names_from=Type,
                       values_from=value)

sep4_lysis=sep4_lysis%>%
  rename(Sample=Time)%>%
  mutate_at("Sample",str_replace,"1","D4_1")%>%
  mutate_at("Sample",str_replace,"3","D4_3")%>%
  mutate_at("Sample",str_replace,"5","D4_5")%>%
  rename(Lysis=R4,
         Phage=P4,
         Host=H4)%>%
  select(Lysis)
#Combine lysis into lysogeny for mat 4
sep4_ll=cbind(sep4_lysogen,sep4_lysis)

#Now combine into master for propr
rna_lys_lys=rbind(sep1_ll,sep2_ll,sep3_ll,sep4_ll)
rna_lys_lys=rna_lys_lys%>%
  select(-Sample)
rna_lys_lys=as.matrix(rna_lys_lys)
#Make proper count matrix for propr
rownames(rna_lys_lys)=rna_lys_lys$Sample
#Now do proportionality anaylsis
library(propr)
propr_lys_lys=propr(counts=rna_lys_lys,
                    metric=c("rho"),
                    ivar="clr",
                    p=100)

snapshot(propr_lys_lys)


###I think because counts are so small (as propr warned) they become too difficult to control for with a clr transform, especially integrating across sample types (RNA and DNA). Lets do this with classic correlation using TPM calculated separately for RNA (all phage genes = 1e6 and the TPM from DNA (phage + host = 1e6)
#Mat 1
sep1=smerged1%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11","R1_2","R1_4"))

#Make new df
llrpk1=sep1
#Put length in KB
llrpk1$Length=llrpk1$Length/1000
#Create RPK for Mat 1
llrpk1$R1_1=llrpk1$R1_1/llrpk1$Length
llrpk1$R1_3=llrpk1$R1_3/llrpk1$Length
llrpk1$R1_5=llrpk1$R1_5/llrpk1$Length
#create scaling factors
#Mat 1
llsum11=(sum(llrpk1$R1_1)/1000000)
llsum13=(sum(llrpk1$R1_3)/1000000)
llsum15=(sum(llrpk1$R1_5)/1000000)
#Create TPM for Mat 1
llrpk1$R1_1=llrpk1$R1_1/llsum11
llrpk1$R1_3=llrpk1$R1_3/llsum13
llrpk1$R1_5=llrpk1$R1_5/llsum15
#Does it make sense
sum(llrpk1$R1_1)
sum(llrpk1$R1_3)
sum(llrpk1$R1_5)

#Bring it back to a sensical name now
D1_ll_tpm=llrpk1

#Extract lysogeny hallmark
D1_lysogeny_tpm_subset=D1_ll_tpm[str_detect(D1_ll_tpm$viral_hit,
                             "integrase | ParA family | ParB | recombinase"),]
#This is not final TPM for phage and host, just merging to get hosts
D1_lysogeny_tpm_subset=merge(D1_lysogeny_tpm_subset,D1_count_tpm,by="Phage")
D1_lysogeny_tpm_subset=D1_lysogeny_tpm_subset%>%
  select(c("R1_1","R1_3","R1_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D1_lysogeny_tpm_subset=pivot_longer(D1_lysogeny_tpm_subset,
                                    cols=2:4,
                                    names_to=c("Type","Time"),
                                    names_sep="_")
D1_lysogeny_tpm_subset=D1_lysogeny_tpm_subset%>%
  rename(Lysogeny=value)%>%
  select(-Type)

#Extract lysis hallmark
D1_lysis_tpm_subset=D1_ll_tpm[str_detect(D1_ll_tpm$viral_hit,
                                            "terminase"),]
D1_lysis_tpm_subset=merge(D1_lysis_tpm_subset,D1_count_tpm,by="Phage")
D1_lysis_tpm_subset=D1_lysis_tpm_subset%>%
  select(c("R1_1","R1_3","R1_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D1_lysis_tpm_subset=pivot_longer(D1_lysis_tpm_subset,
                                    cols=2:4,
                                    names_to=c("Type","Time"),
                                    names_sep="_")
D1_lysis_tpm_subset=D1_lysis_tpm_subset%>%
  rename(Lysis=value)%>%
  select(-Type)

D1_lys_lys_tpm=merge(D1_lysis_tpm_subset,D1_lysogeny_tpm_subset,by=c("Phylum","Time"),all = T)
D1_lys_lys_tpm=D1_lys_lys_tpm%>%
  mutate(Mat="D1")
#End

#Mat 2
sep2=smerged2%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11","R2_2","R2_4"))

#Make new df
llrpk2=sep2
#Put length in KB
llrpk2$Length=llrpk2$Length/1000
#Create RPK for Mat 1
llrpk2$R2_1=llrpk2$R2_1/llrpk2$Length
llrpk2$R2_3=llrpk2$R2_3/llrpk2$Length
llrpk2$R2_5=llrpk2$R2_5/llrpk2$Length
#create scaling factors
#Mat 1
llsum21=(sum(llrpk2$R2_1)/1000000)
llsum23=(sum(llrpk2$R2_3)/1000000)
llsum25=(sum(llrpk2$R2_5)/1000000)
#Create TPM for Mat 1
llrpk2$R2_1=llrpk2$R2_1/llsum21
llrpk2$R2_3=llrpk2$R2_3/llsum23
llrpk2$R2_5=llrpk2$R2_5/llsum25
#Does it make sense
sum(llrpk2$R2_1)
sum(llrpk2$R2_3)
sum(llrpk2$R2_5)

#Bring it back to a sensical name now
D2_ll_tpm=llrpk2

#Extract lysogeny hallmark
D2_lysogeny_tpm_subset=D2_ll_tpm[str_detect(D2_ll_tpm$viral_hit,
                                            "integrase | ParA family | ParB | recombinase"),]
#This is not final TPM for phage and host, just merging to get hosts
D2_lysogeny_tpm_subset=merge(D2_lysogeny_tpm_subset,D2_count_tpm,by="Phage")
D2_lysogeny_tpm_subset=D2_lysogeny_tpm_subset%>%
  select(c("R2_1","R2_3","R2_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D2_lysogeny_tpm_subset=pivot_longer(D2_lysogeny_tpm_subset,
                                    cols=2:4,
                                    names_to=c("Type","Time"),
                                    names_sep="_")
D2_lysogeny_tpm_subset=D2_lysogeny_tpm_subset%>%
  rename(Lysogeny=value)%>%
  select(-Type)

#Extract lysis hallmark
D2_lysis_tpm_subset=D2_ll_tpm[str_detect(D2_ll_tpm$viral_hit,
                                         "terminase"),]
D2_lysis_tpm_subset=merge(D2_lysis_tpm_subset,D2_count_tpm,by="Phage")
D2_lysis_tpm_subset=D2_lysis_tpm_subset%>%
  select(c("R2_1","R2_3","R2_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D2_lysis_tpm_subset=pivot_longer(D2_lysis_tpm_subset,
                                 cols=2:4,
                                 names_to=c("Type","Time"),
                                 names_sep="_")
D2_lysis_tpm_subset=D2_lysis_tpm_subset%>%
  rename(Lysis=value)%>%
  select(-Type)

D2_lys_lys_tpm=merge(D2_lysis_tpm_subset,D2_lysogeny_tpm_subset,by=c("Phylum","Time"),all = T)
D2_lys_lys_tpm=D2_lys_lys_tpm%>%
  mutate(Mat="D2")
#End


#Mat 3
sep3=smerged3%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11","R3_2","R3_4"))

#Make new df
llrpk3=sep3
#Put length in KB
llrpk3$Length=llrpk3$Length/1000
#Create RPK for Mat 1
llrpk3$R3_1=llrpk3$R3_1/llrpk3$Length
llrpk3$R3_3=llrpk3$R3_3/llrpk3$Length
llrpk3$R3_5=llrpk3$R3_5/llrpk3$Length
#create scaling factors
#Mat 1
llsum31=(sum(llrpk3$R3_1)/1000000)
llsum33=(sum(llrpk3$R3_3)/1000000)
llsum35=(sum(llrpk3$R3_5)/1000000)
#Create TPM for Mat 1
llrpk3$R3_1=llrpk3$R3_1/llsum31
llrpk3$R3_3=llrpk3$R3_3/llsum33
llrpk3$R3_5=llrpk3$R3_5/llsum35
#Does it make sense
sum(llrpk3$R3_1)
sum(llrpk3$R3_3)
sum(llrpk3$R3_5)

#Bring it back to a sensical name now
D3_ll_tpm=llrpk3

#Extract lysogeny hallmark
D3_lysogeny_tpm_subset=D3_ll_tpm[str_detect(D3_ll_tpm$viral_hit,
                                            "integrase | ParA family | ParB | recombinase"),]
#This is not final TPM for phage and host, just merging to get hosts
D3_lysogeny_tpm_subset=merge(D3_lysogeny_tpm_subset,D3_count_tpm,by="Phage")
D3_lysogeny_tpm_subset=D3_lysogeny_tpm_subset%>%
  select(c("R3_1","R3_3","R3_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D3_lysogeny_tpm_subset=pivot_longer(D3_lysogeny_tpm_subset,
                                    cols=2:4,
                                    names_to=c("Type","Time"),
                                    names_sep="_")
D3_lysogeny_tpm_subset=D3_lysogeny_tpm_subset%>%
  rename(Lysogeny=value)%>%
  select(-Type)

#Extract lysis hallmark
D3_lysis_tpm_subset=D3_ll_tpm[str_detect(D3_ll_tpm$viral_hit,
                                         "terminase"),]
D3_lysis_tpm_subset=merge(D3_lysis_tpm_subset,D3_count_tpm,by="Phage")
D3_lysis_tpm_subset=D3_lysis_tpm_subset%>%
  select(c("R3_1","R3_3","R3_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D3_lysis_tpm_subset=pivot_longer(D3_lysis_tpm_subset,
                                 cols=2:4,
                                 names_to=c("Type","Time"),
                                 names_sep="_")
D3_lysis_tpm_subset=D3_lysis_tpm_subset%>%
  rename(Lysis=value)%>%
  select(-Type)

D3_lys_lys_tpm=merge(D3_lysis_tpm_subset,D3_lysogeny_tpm_subset,by=c("Phylum","Time"),all = T)
D3_lys_lys_tpm=D3_lys_lys_tpm%>%
  mutate(Mat="D3")
#End


#Mat 4
sep4=smerged4%>%
  separate(col=contig,into=c("1","2","Phage","3","4","5","6","7","8","9","10","11"),sep="_")%>%
  select(-c("1","2","3","4","5","6","7","8","9","10","11","R4_2","R4_4"))

#Make new df
llrpk4=sep4
#Put length in KB
llrpk4$Length=llrpk4$Length/1000
#Create RPK for Mat 1
llrpk4$R4_1=llrpk4$R4_1/llrpk4$Length
llrpk4$R4_3=llrpk4$R4_3/llrpk4$Length
llrpk4$R4_5=llrpk4$R4_5/llrpk4$Length
#create scaling factors
#Mat 1
llsum41=(sum(llrpk4$R4_1)/1000000)
llsum43=(sum(llrpk4$R4_3)/1000000)
llsum45=(sum(llrpk4$R4_5)/1000000)
#Create TPM for Mat 1
llrpk4$R4_1=llrpk4$R4_1/llsum41
llrpk4$R4_3=llrpk4$R4_3/llsum43
llrpk4$R4_5=llrpk4$R4_5/llsum45
#Does it make sense
sum(llrpk4$R4_1)
sum(llrpk4$R4_3)
sum(llrpk4$R4_5)

#Bring it back to a sensical name now
D4_ll_tpm=llrpk4

#Extract lysogeny hallmark
D4_lysogeny_tpm_subset=D4_ll_tpm[str_detect(D4_ll_tpm$viral_hit,
                                            "integrase | ParA family | ParB | recombinase"),]
#This is not final TPM for phage and host, just merging to get hosts
D4_lysogeny_tpm_subset=merge(D4_lysogeny_tpm_subset,D4_count_tpm,by="Phage")
D4_lysogeny_tpm_subset=D4_lysogeny_tpm_subset%>%
  select(c("R4_1","R4_3","R4_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D4_lysogeny_tpm_subset=pivot_longer(D4_lysogeny_tpm_subset,
                                    cols=2:4,
                                    names_to=c("Type","Time"),
                                    names_sep="_")
D4_lysogeny_tpm_subset=D4_lysogeny_tpm_subset%>%
  rename(Lysogeny=value)%>%
  select(-Type)

#Extract lysis hallmark
D4_lysis_tpm_subset=D4_ll_tpm[str_detect(D4_ll_tpm$viral_hit,
                                         "terminase"),]
D4_lysis_tpm_subset=merge(D4_lysis_tpm_subset,D4_count_tpm,by="Phage")
D4_lysis_tpm_subset=D4_lysis_tpm_subset%>%
  select(c("R4_1","R4_3","R4_5","Phylum"))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
D4_lysis_tpm_subset=pivot_longer(D4_lysis_tpm_subset,
                                 cols=2:4,
                                 names_to=c("Type","Time"),
                                 names_sep="_")
D4_lysis_tpm_subset=D4_lysis_tpm_subset%>%
  rename(Lysis=value)%>%
  select(-Type)

D4_lys_lys_tpm=merge(D4_lysis_tpm_subset,D4_lysogeny_tpm_subset,by=c("Phylum","Time"),all = T)
D4_lys_lys_tpm=D4_lys_lys_tpm%>%
  mutate(Mat="D4")
#End

#Merge into one master frame
lys_lys_tpm_1=rbind(D1_lys_lys_tpm,D2_lys_lys_tpm,D3_lys_lys_tpm,D4_lys_lys_tpm)
lys_lys_tpm_1=lys_lys_tpm_1%>%
  mutate_at("Phylum",str_replace,"p__","")%>%
  mutate_at("Phylum",str_replace,"_F","")
#Merge with p and h counts
lys_lys_tpm=merge(lys_lys_tpm_1,Phage_host_count_merge,by=c("Phylum","Time","Mat"))
#Create abundance corrected lysis and lysogeny, and log values of each
lys_lys_tpm=lys_lys_tpm%>%
  mutate(abundLysis=Lysis/P)%>%
  mutate(abundLysogeny=Lysogeny/P)%>%
  mutate(logLysis=log10(Lysis+1))%>%
  mutate(logLysogeny=log10(Lysogeny+1))%>%
  mutate(logAbundLysis=log10(abundLysis+1))%>%
  mutate(logAbundLysogeny=log10(abundLysogeny+1))

#Lets get to modelling
#Try Dan model
conflict_prefer("lmer", "lme4")
count_glm_lysogeny=lmer(data=lys_lys_tpm,logLysogeny~logVMR+(1+logVMR|Phylum),na.action = na.omit)
summary(count_glm_lysogeny)
anova(count_glm_lysogeny)
fixef(count_glm_lysogeny,condVars=T)
ranef(count_glm_lysogeny,condVar=T)
coef(count_glm_lysogeny)
confint(count_glm_lysogeny)

#Test assumptions
library("DHARMa")
check_vmr_count_model_lysogeny <- simulateResiduals(fittedModel = count_glm_lysogeny, n = 999)
plot(check_vmr_count_model_lysogeny)
par(mfrow=c(2,2))
plot(count_glm_lysogeny)

#Create dataframe of range for passing to geom line
line_lysogeny_df=data.frame(VMR=c(min(lys_lys_tpm$logVMR),
                            max(lys_lys_tpm$logVMR)),
                   lysogeny=c(fixef(count_glm_lysogeny)[1] +
                              fixef(count_glm_lysogeny)[2]*
                              min(lys_lys_tpm$logVMR),
                            fixef(count_glm_lysogeny)[1] +
                              fixef(count_glm_lysogeny)[2]*
                              max(lys_lys_tpm$logVMR)))
#plot lysogeny
ggplot(lys_lys_tpm, aes(y=logLysogeny,x=logVMR))+
  geom_point(alpha=0.8, aes(color=Phylum))+
  geom_line(data=line_lysogeny_df,
            aes(x=VMR,y=lysogeny),
            linetype="solid",
            color="#333333",
            size=1)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Gene Expression ("~log[10]~TPM~')'),
       x=expression("VMR ("~log[10]~')'))+
  scale_color_viridis_d(direction=-1)
dev.print(png,fil="vmr_overall_phyla.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)


#Lysis model
conflict_prefer("lmer", "lme4")
count_glm_lysis=lmer(data=lys_lys_tpm,logLysis~logVMR+(1+logVMR|Phylum))
summary(count_glm_lysis)
anova(count_glm_lysis)
fixef(count_glm_lysis,condVars=T)
ranef(count_glm_lysis,condVar=T)
coef(count_glm_lysis)
confint(count_glm_lysis)

#Test assumptions
library("DHARMa")
check_vmr_count_model_lysis <- simulateResiduals(fittedModel = count_glm_lysis, n = 999)
plot(check_vmr_count_model_lysis)
par(mfrow=c(2,2))
plot(count_glm_lysis)

#Nothin is significant, what about just time
#Time and lysogeny
conflict_prefer("lmer", "lme4")
time_glm_lysogeny=glm(data=lys_lys_tpm,logLysogeny~Time,na.action = na.omit)
summary(time_glm_lysogeny)
Anova(time_glm_lysogeny,Type="III")
#Nope
fixef(time_glm_lysogeny,condVars=T)
ranef(time_glm_lysogeny,condVar=T)
coef(time_glm_lysogeny)
confint(time_glm_lysogeny)
library(multcomp)
summary(glht(time_glm_lysogeny, linfct=mcp(Time="Tukey")))
#Test assumptions
library("DHARMa")
check_vmr_count_model_lysogeny <- simulateResiduals(fittedModel = time_glm_lysogeny, n = 999)
plot(check_vmr_count_model_lysogeny)
par(mfrow=c(2,2))
plot(time_glm_lysogeny)

#plot Lysogeny across time
ggplot(lys_lys_tpm, aes(y=logLysogeny,x=Time))+
  geom_boxplot(aes(fill=Time),
               outlier.shape = NA,
               size=1,
               alpha=0.9)+
  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2"),name="Time")+
  geom_point(aes(fill=Time),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMR ("~log[10]~')'),
       x="Sampling Time Point")

#Nothin is significant, what about just time
#Time and lysis
conflict_prefer("lmer", "lme4")
time_glm_lysis=glm(data=lys_lys_tpm,logLysis~Time,na.action = na.omit)
summary(time_glm_lysis)
Anova(time_glm_lysis,Type="III")
#Nope
#Test assumptions
library("DHARMa")
check_vmr_count_model_lysis <- simulateResiduals(fittedModel = time_glm_lysis, n = 999)
plot(check_vmr_count_model_lysis)
par(mfrow=c(2,2))
plot(time_glm_lysogeny)

ggplot(lys_lys_tpm, aes(y=logLysis,x=Time))+
  geom_boxplot(aes(fill=Time),
               outlier.shape = NA,
               size=1,
               alpha=0.9)+
  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2"),name="Time")+
  geom_point(aes(fill=Time),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMR ("~log[10]~')'),
       x="Sampling Time Point")



#################
#################
#################
#################
#################
#################
#################
####TPM plots for interpretable visualization of expression including both DNA and RNA kept to MAG-resolved level
#read in DNA data
dna_diel_region1=read.xls("/Users/cissell/Desktop/D1_DNA_region_counts.xlsx")
dna_diel_region2=read.xls("/Users/cissell/Desktop/D2_DNA_region_counts.xlsx")
dna_diel_region3=read.xls("/Users/cissell/Desktop/D3_DNA_region_counts.xlsx")
dna_diel_region4=read.xls("/Users/cissell/Desktop/D4_DNA_region_counts.xlsx")
#Merge rna counts with gene ids
merged1=merge(rna_diel1,rna_kegg1, by="contig",all=T)
merged2=merge(rna_diel2,rna_kegg2, by="contig",all=T)
merged3=merge(rna_diel3,rna_kegg3, by="contig",all=T)
merged4=merge(rna_diel4,rna_kegg4, by="contig",all=T)
#Merge in DNA counts
merged1=merge(merged1,dna_diel_region1,all=T)
merged2=merge(merged2,dna_diel_region2,all=T)
merged3=merge(merged3,dna_diel_region3,all=T)
merged4=merge(merged4,dna_diel_region4,all=T)

#Remove rows without KEGG ID
smerged1=merged1 %>%
  na_if("") %>%
  na.omit
smerged2=merged2 %>%
  na_if("") %>%
  na.omit
smerged3=merged3 %>%
  na_if("") %>%
  na.omit
smerged4=merged4 %>%
  na_if("") %>%
  na.omit
#Drop contig
smerged1=subset(smerged1,select=-c(contig))
smerged2=subset(smerged2,select=-c(contig))
smerged3=subset(smerged3,select=-c(contig))
smerged4=subset(smerged4,select=-c(contig))

#Merge to subset
tpmmerged1=merge(energy1,smerged1, by="kegg_id")
tpmmerged2=merge(energy2,smerged2, by="kegg_id")
tpmmerged3=merge(energy3,smerged3, by="kegg_id")
tpmmerged4=merge(energy4,smerged4, by="kegg_id")

#Add back in SQR
#Subset out by kegg
sqr1=smerged1%>%
  filter(kegg_id=="K17218")
sqr2=smerged2%>%
  filter(kegg_id=="K17218")
sqr3=smerged3%>%
  filter(kegg_id=="K17218")
sqr4=smerged4%>%
  filter(kegg_id=="K17218")
#then combine
tpmmerged1=rbind.fill(tpmmerged1,sqr1)
tpmmerged2=rbind.fill(tpmmerged2,sqr2)
tpmmerged3=rbind.fill(tpmmerged3,sqr3)
tpmmerged4=rbind.fill(tpmmerged4,sqr4)

#Make new df
rpk1=tpmmerged1
rpk2=tpmmerged2
rpk3=tpmmerged3
rpk4=tpmmerged4
#Put length in KB
rpk1$Length=rpk1$Length/1000
rpk2$Length=rpk2$Length/1000
rpk3$Length=rpk3$Length/1000
rpk4$Length=rpk4$Length/1000
#Create RPK for Mat 1
rpk1$R1_1=rpk1$R1_1/rpk1$Length
rpk1$R1_2=rpk1$R1_2/rpk1$Length
rpk1$R1_3=rpk1$R1_3/rpk1$Length
rpk1$R1_4=rpk1$R1_4/rpk1$Length
rpk1$R1_5=rpk1$R1_5/rpk1$Length
rpk1$D1_1=rpk1$D1_1/rpk1$Length
rpk1$D1_3=rpk1$D1_3/rpk1$Length
rpk1$D1_5=rpk1$D1_5/rpk1$Length
#Create RPK for Mat 2
rpk2$R2_1=rpk2$R2_1/rpk2$Length
rpk2$R2_2=rpk2$R2_2/rpk2$Length
rpk2$R2_3=rpk2$R2_3/rpk2$Length
rpk2$R2_4=rpk2$R2_4/rpk2$Length
rpk2$R2_5=rpk2$R2_5/rpk2$Length
rpk2$D2_1=rpk2$D2_1/rpk2$Length
rpk2$D2_3=rpk2$D2_3/rpk2$Length
rpk2$D2_5=rpk2$D2_5/rpk2$Length
#Create RPK for Mat 3
rpk3$R3_1=rpk3$R3_1/rpk3$Length
rpk3$R3_2=rpk3$R3_2/rpk3$Length
rpk3$R3_3=rpk3$R3_3/rpk3$Length
rpk3$R3_4=rpk3$R3_4/rpk3$Length
rpk3$R3_5=rpk3$R3_5/rpk3$Length
rpk3$D3_1=rpk3$D3_1/rpk3$Length
rpk3$D3_3=rpk3$D3_3/rpk3$Length
rpk3$D3_5=rpk3$D3_5/rpk3$Length
#Create RPK for Mat 4
rpk4$R4_1=rpk4$R4_1/rpk4$Length
rpk4$R4_2=rpk4$R4_2/rpk4$Length
rpk4$R4_3=rpk4$R4_3/rpk4$Length
rpk4$R4_4=rpk4$R4_4/rpk4$Length
rpk4$R4_5=rpk4$R4_5/rpk4$Length
rpk4$D4_1=rpk4$D4_1/rpk4$Length
rpk4$D4_3=rpk4$D4_3/rpk4$Length
rpk4$D4_5=rpk4$D4_5/rpk4$Length
#create scaling factors
#Mat 1
sum11=(sum(rpk1$R1_1)/1000000)
sum12=(sum(rpk1$R1_2)/1000000)
sum13=(sum(rpk1$R1_3)/1000000)
sum14=(sum(rpk1$R1_4)/1000000)
sum15=(sum(rpk1$R1_5)/1000000)
sumD11=(sum(rpk1$D1_1)/1000000)
sumD13=(sum(rpk1$D1_3)/1000000)
sumD15=(sum(rpk1$D1_5)/1000000)
#Mat 2
sum21=(sum(rpk2$R2_1)/1000000)
sum22=(sum(rpk2$R2_2)/1000000)
sum23=(sum(rpk2$R2_3)/1000000)
sum24=(sum(rpk2$R2_4)/1000000)
sum25=(sum(rpk2$R2_5)/1000000)
sumD21=(sum(rpk2$D2_1)/1000000)
sumD23=(sum(rpk2$D2_3)/1000000)
sumD25=(sum(rpk2$D2_5)/1000000)
#Mat 3
sum31=(sum(rpk3$R3_1)/1000000)
sum32=(sum(rpk3$R3_2)/1000000)
sum33=(sum(rpk3$R3_3)/1000000)
sum34=(sum(rpk3$R3_4)/1000000)
sum35=(sum(rpk3$R3_5)/1000000)
sumD31=(sum(rpk3$D3_1)/1000000)
sumD33=(sum(rpk3$D3_3)/1000000)
sumD35=(sum(rpk3$D3_5)/1000000)
#Mat 4
sum41=(sum(rpk4$R4_1)/1000000)
sum42=(sum(rpk4$R4_2)/1000000)
sum43=(sum(rpk4$R4_3)/1000000)
sum44=(sum(rpk4$R4_4)/1000000)
sum45=(sum(rpk4$R4_5)/1000000)
sumD41=(sum(rpk4$D4_1)/1000000)
sumD43=(sum(rpk4$D4_3)/1000000)
sumD45=(sum(rpk4$D4_5)/1000000)
#Create TPM for Mat 1
rpk1$R1_1=rpk1$R1_1/sum11
rpk1$R1_2=rpk1$R1_2/sum12
rpk1$R1_3=rpk1$R1_3/sum13
rpk1$R1_4=rpk1$R1_4/sum14
rpk1$R1_5=rpk1$R1_5/sum15
rpk1$D1_1=rpk1$D1_1/sumD11
rpk1$D1_3=rpk1$D1_3/sumD13
rpk1$D1_5=rpk1$D1_5/sumD15
#Create TPM for Mat 2
rpk2$R2_1=rpk2$R2_1/sum21
rpk2$R2_2=rpk2$R2_2/sum22
rpk2$R2_3=rpk2$R2_3/sum23
rpk2$R2_4=rpk2$R2_4/sum24
rpk2$R2_5=rpk2$R2_5/sum25
rpk2$D2_1=rpk2$D2_1/sumD21
rpk2$D2_3=rpk2$D2_3/sumD23
rpk2$D2_5=rpk2$D2_5/sumD25
#Create TPM for Mat 3
rpk3$R3_1=rpk3$R3_1/sum31
rpk3$R3_2=rpk3$R3_2/sum32
rpk3$R3_3=rpk3$R3_3/sum33
rpk3$R3_4=rpk3$R3_4/sum34
rpk3$R3_5=rpk3$R3_5/sum35
rpk3$D3_1=rpk3$D3_1/sumD31
rpk3$D3_3=rpk3$D3_3/sumD33
rpk3$D3_5=rpk3$D3_5/sumD35
#Create TPM for Mat 4
rpk4$R4_1=rpk4$R4_1/sum41
rpk4$R4_2=rpk4$R4_2/sum42
rpk4$R4_3=rpk4$R4_3/sum43
rpk4$R4_4=rpk4$R4_4/sum44
rpk4$R4_5=rpk4$R4_5/sum45
rpk4$D4_1=rpk4$D4_1/sumD41
rpk4$D4_3=rpk4$D4_3/sumD43
rpk4$D4_5=rpk4$D4_5/sumD45

#Check for sense
sum(rpk1$R1_1) #S'all good
sum(rpk2$D2_5) #S'all good
sum(rpk3$D3_3) #S'all good

#Merge across KEGG ID and Class
conflict_prefer("rename", "dplyr")
tpm1=rpk1 %>%
  group_by(kegg_id,Mag) %>%
  summarise(across(c(D1_1,D1_3,D1_5,R1_1,R1_2,R1_3,R1_4,R1_5), sum))%>%
  ungroup()%>%
  rename(MAG=Mag)%>%
  group_by(kegg_id,MAG)%>%
  mutate(D=mean(c(D1_1,D1_3,D1_5)))%>%
  mutate(R=mean(c(R1_1,R1_2,R1_3,R1_4,R1_5)))%>%
  select(-c(R1_1,R1_2,R1_3,R1_4,R1_5,D1_1,D1_3,D1_5))

tpm2=rpk2 %>%
  group_by(kegg_id,MAG) %>%
  summarise(across(c(D2_1,D2_3,D2_5,R2_1,R2_2,R2_3,R2_4,R2_5), sum))%>%
  ungroup()%>%
  group_by(kegg_id,MAG)%>%
  mutate(D=mean(c(D2_1,D2_3,D2_5)))%>%
  mutate(R=mean(c(R2_1,R2_2,R2_3,R2_4,R2_5)))%>%
  select(-c(R2_1,R2_2,R2_3,R2_4,R2_5,D2_1,D2_3,D2_5))

tpm3=rpk3 %>%
  group_by(kegg_id,MAG) %>%
  summarise(across(c(D3_1,D3_3,D3_5,R3_1,R3_2,R3_3,R3_4,R3_5), sum))%>%
  ungroup()%>%
  group_by(kegg_id,MAG)%>%
  mutate(D=mean(c(D3_1,D3_3,D3_5)))%>%
  mutate(R=mean(c(R3_1,R3_2,R3_3,R3_4,R3_5)))%>%
  select(-c(R3_1,R3_2,R3_3,R3_4,R3_5,D3_1,D3_3,D3_5))

tpm4=rpk4 %>%
  group_by(kegg_id,MAG) %>%
  summarise(across(c(D4_1,D4_3,D4_5,R4_1,R4_2,R4_3,R4_4,R4_5), sum))%>%
  ungroup()%>%
  group_by(kegg_id,MAG)%>%
  mutate(D=mean(c(D4_1,D4_3,D4_5)))%>%
  mutate(R=mean(c(R4_1,R4_2,R4_3,R4_4,R4_5)))%>%
  select(-c(R4_1,R4_2,R4_3,R4_4,R4_5,D4_1,D4_3,D4_5))

#Merge to get descriptions again
tpmmerged1=merge(energy1,tpm1, by="kegg_id",all.y=T)
tpmmerged2=merge(energy2,tpm2, by="kegg_id",all.y=T)
tpmmerged3=merge(energy3,tpm3, by="kegg_id",all.y=T)
tpmmerged4=merge(energy4,tpm4, by="kegg_id",all.y=T)

#Remove redundancy
tpmmerged1=distinct(tpmmerged1,kegg_id,MAG,.keep_all = T)
tpmmerged2=distinct(tpmmerged2,kegg_id,MAG,.keep_all = T)
tpmmerged3=distinct(tpmmerged3,kegg_id,MAG,.keep_all = T)
tpmmerged4=distinct(tpmmerged4,kegg_id,MAG,.keep_all = T)

#Merge across gene ids
tpm_merge1=merge(y=tpmmerged1,x=tpmmerged2,by=c("kegg_id","MAG","module","gene_description","D","R"),all=T)
tpm_merge2=merge(tpmmerged3,tpmmerged4,by=c("kegg_id","MAG","module","gene_description","D","R"),all=T)
tpm_merge=merge(y=tpm_merge1,x=tpm_merge2,by=c("kegg_id","MAG","module","gene_description","D","R"),all=T)

#Replace NA with 0
tpm_merge$D[is.na(tpm_merge$D)] <- 0
tpm_merge$R[is.na(tpm_merge$R)] <- 0

#Subset out important genes (Carbon, Sulfur, Nitrogen, Phosphorous)
tpm_merge_subs=tpm_merge[str_detect(tpm_merge$kegg_id,"K02703|K02706|K02705|K02689|K02694|K02690|K02691|K02692|K02274|K02275|K02276|K02277|K17218|K13811|K00958|K00956|K00957|K00860|K00390|K00380|K00381|K00392|K00394|K00395|K11180|K11181|K17222|K17223|K17226|K17227|K17224|K17225|K22622|K02588|K10946|K10944|K10945|K10535|K00371|K00370|K00374|K02567|K02568|K00368|K15864|K04561|K02305|K00376|K00362|K00363|K03385|K15876|K00367|K10534|K00372|K00360|K17877|K00366|K20932|K20933|K20934|K20935|K04758"),]
conflict_prefer("filter", "dplyr")
#tpm_long_subs=tpm_long_subs%>%
#  filter(!(Class=="Bacteroidia" & module=="Photosystem II"))
#Plot
#Set factor order of Kegg ID for plotting
tpm_merge_subs$kegg_id=factor(tpm_merge_subs$kegg_id,levels=c("K02703","K02706","K02705","K02689","K02694","K02690","K02691","K02692","K02274","K02275","K02276","K02277","K17218","K13811","K00958","K00956","K00957","K00860","K00390","K00380","K00381","K00392","K00394","K00395","K11180","K11181","K17222","K17223","K17226","K17227","K17224","K17225","K22622","K02588","K10946","K10944","K10945","K10535","K00370","K00371","K00374","K02567","K02568","K00368K15864","K04561","K02305","K00376","K00362","K00363","K03385","K15876","K00367","K10534","K00372","K00360","K17877","K00366","K20932","K20933","K20934","K20935","K04758"))
tpm_merge_subs$MAG=factor(tpm_merge_subs$MAG,levels=c("D2_MAG_00001","D1_MAG_00001","D4_MAG_00001","D3_MAG_00019","D1_MAG_00047","D2_MAG_00034","D1_MAG_00004","D2_MAG_00002","D1_MAG_00049","D2_MAG_00047","D4_MAG_00019","D1_MAG_00031","D1_MAG_00022","D4_MAG_00008","D3_MAG_00040","D3_MAG_00042","D2_MAG_00028","D4_MAG_00013","D1_MAG_00025","D4_MAG_00045","D4_MAG_00029","D1_MAG_00037","D3_MAG_00007","D1_MAG_00035","D4_MAG_00025","D2_MAG_00007","D3_MAG_00031","D3_MAG_00078","D3_MAG_00041","D1_MAG_00016","D2_MAG_00062","D2_MAG_00061","D3_MAG_00036","D3_MAG_00015","D1_MAG_00034","D2_MAG_00022","D4_MAG_00052","D3_MAG_00016","D1_MAG_00044","D3_MAG_00065","D2_MAG_00027","D4_MAG_00040","D2_MAG_00056","D4_MAG_00023","D1_MAG_00008","D4_MAG_00006","D3_MAG_00012","D2_MAG_00009","D3_MAG_00010","D2_MAG_00042","D2_MAG_00055","D2_MAG_00064","D2_MAG_00060","D2_MAG_00066","D2_MAG_00065","D3_MAG_00075","D4_MAG_00048","D3_MAG_00069","D2_MAG_00046","D2_MAG_00024","D1_MAG_00012","D4_MAG_00011","D3_MAG_00072","D3_MAG_00028","D4_MAG_00058","D3_MAG_00052","D1_MAG_00030","D2_MAG_00051","D3_MAG_00035","D4_MAG_00041","D2_MAG_00015","D3_MAG_00033","D1_MAG_00011","D3_MAG_00018","D3_MAG_00001","D2_MAG_00012","D1_MAG_00017","D3_MAG_00059","D2_MAG_00053","D4_MAG_00014","D3_MAG_00008","D2_MAG_00033","D1_MAG_00023","D1_MAG_00018","D3_MAG_00002","D2_MAG_00017","D3_MAG_00074","D1_MAG_00050","D4_MAG_00039","D3_MAG_00006","D4_MAG_00051","D3_MAG_00071","D2_MAG_00023","D3_MAG_00073","D2_MAG_00011","D1_MAG_00036","D2_MAG_00010","D3_MAG_00053","D1_MAG_00027","D3_MAG_00034","D1_MAG_00045","D2_MAG_00045","D1_MAG_00041","D4_MAG_00033","D3_MAG_00058","D4_MAG_00055","D3_MAG_00054","D2_MAG_00040","D1_MAG_00053","D3_MAG_00077","D2_MAG_00031","D4_MAG_00044","D1_MAG_00026","D2_MAG_00048","D4_MAG_00027","D3_MAG_00045","D4_MAG_00004","D3_MAG_00024","D1_MAG_00043","D2_MAG_00038","D4_MAG_00017","D1_MAG_00013","D3_MAG_00009","D4_MAG_00012","D2_MAG_00026","D4_MAG_00016","D1_MAG_00029","D3_MAG_00061","D2_MAG_00052","D1_MAG_00048","D2_MAG_00013","D1_MAG_00010","D4_MAG_00015","D3_MAG_00023","D4_MAG_00030","D3_MAG_00049","D4_MAG_00049","D3_MAG_00068","D3_MAG_00051","D2_MAG_00043","D4_MAG_00034","D4_MAG_00026","D3_MAG_00021","D1_MAG_00020","D3_MAG_00020","D2_MAG_00063","D4_MAG_00035","D4_MAG_00043","D3_MAG_00030","D2_MAG_00050","D1_MAG_00038","D4_MAG_00007","D3_MAG_00013","D2_MAG_00019","D1_MAG_00028","D2_MAG_00049","D4_MAG_00038","D1_MAG_00056","D3_MAG_00060","D1_MAG_00024","D3_MAG_00055","D1_MAG_00006","D3_MAG_00047","D1_MAG_00021","D3_MAG_00038","D4_MAG_00009","D2_MAG_00039","D1_MAG_00019","D4_MAG_00042","D1_MAG_00040","D4_MAG_00036","D3_MAG_00057","D2_MAG_00044","D3_MAG_00066","D2_MAG_00041","D1_MAG_00054","D3_MAG_00043","D4_MAG_00031","D3_MAG_00046","D2_MAG_00006","D2_MAG_00016","D2_MAG_00057","D3_MAG_00048","D4_MAG_00032","D3_MAG_00039","D2_MAG_00020","D4_MAG_00020","D1_MAG_00014","D4_MAG_00050","D1_MAG_00055","D3_MAG_00064","D1_MAG_00039","D4_MAG_00056","D3_MAG_00076","D1_MAG_00052","D3_MAG_00063","D3_MAG_00056","D2_MAG_00059","D3_MAG_00067","D3_MAG_00022","D2_MAG_00029","D1_MAG_00015","D1_MAG_00003","D2_MAG_00005","D4_MAG_00002","D3_MAG_00004","D2_MAG_00025","D1_MAG_00002","D3_MAG_00003","D2_MAG_00003","D4_MAG_00010","D4_MAG_00021","D4_MAG_00028","D3_MAG_00011","D2_MAG_00037","D4_MAG_00024","D3_MAG_00017","D1_MAG_00051","D3_MAG_00044","D1_MAG_00032","D2_MAG_00021","D4_MAG_00005","D3_MAG_00014","D3_MAG_00062","D3_MAG_00027","D4_MAG_00037","D3_MAG_00026","D2_MAG_00014","D1_MAG_00005","D4_MAG_00022","D2_MAG_00030","D2_MAG_00018","D1_MAG_00007","D4_MAG_00059","D2_MAG_00032","D3_MAG_00050","D1_MAG_00046","D4_MAG_00053","D2_MAG_00067","D2_MAG_00008","D3_MAG_00029","D1_MAG_00057","D1_MAG_00009","D4_MAG_00003","D3_MAG_00005","D2_MAG_00004","D3_MAG_00032","D2_MAG_00058","D2_MAG_00054","D3_MAG_00070","D3_MAG_00037","D4_MAG_00046","D4_MAG_00047","D4_MAG_00057","D1_MAG_00033","D4_MAG_00018","D2_MAG_00035","D3_MAG_00025","D2_MAG_00036"))


#Filter out errant classification
tpm_merge_subs=tpm_merge_subs%>%
  filter(!(kegg_id=="K02691" &
             MAG=="D3_MAG_00015"))%>%
  filter(!(kegg_id=="K02691" &
             MAG=="D3_MAG_00036"))%>%
  filter(!(kegg_id=="K02691" &
             MAG=="D4_MAG_00057"))
tpm_merge_subs=tpm_merge_subs%>%
  filter(!(kegg_id=="K02706" &
             MAG=="D4_MAG_00037"))%>%
  filter(!(kegg_id=="K02705" & MAG=="D4_MAG_00037"))
str(tpm_merge_subs)

tpm_merge_subs=tpm_merge_subs[!is.na(tpm_merge_subs$kegg_id),]

#Add in sample ID for finding novel functions within mat
tpm_merge_subs_novel=tpm_merge_subs %>%
  mutate(Mat=case_when(
    str_detect(MAG,"D1") ~ "1",
    str_detect(MAG,"D2") ~ "2",
    str_detect(MAG,"D3") ~ "3",
    str_detect(MAG,"D4") ~ "4"))

diel_rna_heat=ggplot(data=tpm_merge_subs,aes(x=kegg_id,
                                            y=MAG,
                                            fill=log10(1+R)))+
  geom_tile(stat="identity",width=0.5,height=0.95)+
  theme_classic()+
  scale_fill_gradient(low="white",high="#3d5a80")+
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1))+
  theme(plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"))+
  labs(fill="RNA(TPM)")
diel_rna_heat
ggsave(diel_rna_heat,filename="/Users/cissell/Desktop/rna.svg",bg="transparent")

diel_dna_heat=ggplot(data=tpm_merge_subs,aes(x=kegg_id,
                                             y=MAG,
                                             fill=log10(1+D)))+
  geom_tile(stat="identity",width=0.5,height=0.95)+
  theme_classic()+
  scale_fill_gradient(low="white",high="#ee6c4d")+
  theme(axis.text.x=element_text(angle=90,
                                 vjust=1,
                                 hjust=1))+
  theme(plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"))+
  labs(fill="DNA(TPM)")
diel_dna_heat
ggsave(diel_dna_heat,filename="/Users/cissell/Desktop/dna.svg",bg="transparent")


###############################
