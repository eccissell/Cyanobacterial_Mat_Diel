library(dplyr)
library(ggplot2)
theme_set(theme_classic())
library(viridis)
library(tidyverse)
library(gdata)
library(reshape2)
library(phyloseq)
library(vegan)
library(ggpubr)
library(conflicted)
library(codyn)
library(binom)
#Diel Sequencing DNA ORDER data from Kaiju (reads)
##########################
##Lets try using the raw reads with a centered-log-ratio transformation and compare
#####DIEL SEQUENCING DATA
##Read in data
diel=read.xls('/Users/cissell/Desktop/kaiju_summary.xlsx')

##Wrangle the data
otu=pivot_wider(diel,id_cols=taxon_name,names_from = file, values_from = reads)
#Replace NA with 0
otu[is.na(otu)] <- 0

#Remove Unclassified
otu=otu[-c(64),]
otu=otu[-c(62),]
otu=otu[-c(62),]

#Make OTU table for phyloseq
#Rownames have to be ID
rownames(otu)<- otu$taxon_name


#Check counts
count <- colSums(otu[, c(2:ncol(otu))])
count

#Create taxonomy table
taxa<- otu %>% 
  select(taxon_name) %>% 
  separate(taxon_name, c("Domain", "Phylum", "Class", "Order"),
           ";")
str(taxa)

#Make factors not strings
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
str(taxa)
#Add first column of OTU table so OTU and Taxa tables match
taxa<- cbind(otu$taxon_name, taxa)
#Rename first column
colnames(taxa)[1]<- "taxon_name"
str(taxa)
#OTU IDs have to be rownames, just like OTU table
rownames(taxa)<- taxa$taxon_name

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
otu<- subset(otu, select=-taxon_name)
rownames(otu)<- taxa$taxon_name
taxa<- taxa %>% 
  select(-taxon_name)

#Create metadata table
meta=read.xls('/Users/cissell/Desktop/metadata_table_diel.xlsx')
str(meta)
rownames(meta)<- meta$SampleID
meta<- meta %>% 
  select(-SampleID)
#Make proper levels
meta$Mat<- factor(meta$Mat, levels = c("D1", "D2", "D3","D4"))
meta$Age<- factor(meta$Age, levels = c("1", "3", "5"))

#Create matrices before creating phyloseq object
otu_mat<- as.matrix(otu)
otu_mat=otu_mat
tax_mat<- tax_table(as.matrix(taxa))

#Transform to phyloseq objects
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples<- sample_data(meta)

#And unite them into one object
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)
#Check and see if it makes sense
sample_sums(phylo_object)
sample_names(phylo_object)
rank_names(phylo_object) 
sample_variables(phylo_object) 
otu_table(phylo_object)[1:3, 1:2]
taxa_names(phylo_object)[1:5]
tax_table(phylo_object)
taxa_names(phylo_object)

#transform the data using a centered-log-ratio transformation
library(microbiome)
phylo_clr=microbiome::transform(phylo_object, transform='clr',target="sample")
#otu_table(phylo_clr)[otu_table(phylo_clr) < 0.0] <- 0.0

#NMDS Analysis
phylo_clr_nmds<- ordinate(phylo_clr, method = "NMDS", distance="euclidean")
phylo_clr_nmds
#Shepards test/goodness of fit NMDS
goodness(phylo_clr_nmds)
stressplot(phylo_clr_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)
dev.print(png,fil="/Users/cissell/Desktop/stress.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

#Create score matrix
data.scores=as.data.frame(scores(phylo_clr_nmds))
agev=c("1","3","5","1","3","5","1","3","5","1","3","5")
matv=c("1","1","1","2","2","2","3","3","3","4","4","4")
data.scores$Age=agev
data.scores$Mat=matv
#Make hull vectors and hull data frame for plotting polygons
age.1 = data.scores[data.scores$Age == "1",][chull(data.scores[data.scores$Age=="1",c("NMDS1","NMDS2")]),]
age.2=data.scores[data.scores$Age == "3",][chull(data.scores[data.scores$Age=="3",c("NMDS1","NMDS2")]),]
age.3=data.scores[data.scores$Age == "5",][chull(data.scores[data.scores$Age=="5",c("NMDS1","NMDS2")]),]
hull.data=rbind(age.1,age.2,age.3)
hull.data
#Plot NMDS
nmds=plot_ordination(phylo_clr, phylo_clr_nmds,
                    color = "Age", shape = "Mat") +
  theme_classic() +geom_point(size=5.9)+
  scale_color_manual(values=c("#90a295","#455765","#cd9fb2"))+
  scale_shape_manual(values=c(19,17,15,18))+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1))+
  guides(shape=guide_legend(override.aes=list(size=5)))
  
#nmds = nmds +
#  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2))

nmds


ggplot()+
  geom_text(data=species.scores,aes)

#Try plotting from raw scores
nmds2=ggplot() +
  geom_polygon(data=hull.data,
               aes(x=NMDS1,y=NMDS2,fill=Age,group=Age),
               alpha=0.4)+
  geom_polygon(data=data.scores,
               aes(x=NMDS1,y=NMDS2, group=Mat),
               size=1,
               linetype=3, 
               fill=NA,
               colour="#333333")+
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, shape=Mat,colour=Age),size=6) +
  scale_color_manual(values=c("#90a295","#455765","#cd9fb2"),name="Time") +
  scale_fill_manual(values=c("1"="#90a295","3"="#455765","5"="#cd9fb2")) +
  scale_shape_manual(values=c(19,17,15,18))+
  coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=1.3,y=-1.5, label="Stress = 0.125", color="black",size=5)+
  annotate(geom="text",x=0,y=-1, label="3", color="#333333",size=5)+
  annotate(geom="text",x=0.78,y=-0.42, label="1", color="#333333",size=5)+
  annotate(geom="text",x=-0.8,y=0.2, label="4", color="#333333",size=5)+
  annotate(geom="text",x=0.3,y=1.25, label="2", color="#333333",size=5)
nmds2

dev.print(png,fil="nmds2.png",type="quartz",antialias="default",width=6.5,
          height=5.5,units="in",res=1300)


#NMDS on raw dataframe for species overlays
#MAke dataframe
nmds_df=psmelt(phylo_clr)
View(nmds_df)
nmds_df=nmds_df%>%
  select(-OTU,-Domain,-Phylum,-Class)%>%
  spread(nmds_df,key=Order,value=Abundance)

dat=nmds_df[,4:77]
dat=as.data.frame(lapply(dat,unlist))
str(dat)
grp=as.data.frame(nmds_df[,2:3],header=FALSE)
new_nmds=metaMDS(dat,distance="euclidean")
new_nmds

stressplot(new_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)

order.fit=envfit(new_nmds,dat,permutations=999)
plot(new_nmds,type="t")
data.score=as.data.frame(scores(new_nmds))
data.score$site=rownames(data.score)
data.score$grp=grp

species.scores=envfit(new_nmds,dat,permutations=999)
head(species.scores)
species.scores2=as.data.frame(scores(species.scores,display="vectors"))
species.scores2=cbind(species.scores2, Orders=rownames(species.scores2))
species.scores2=cbind(species.scores2,pval=species.scores$vectors$pvals)
head(species.scores2)
sig.sp.sc=subset(species.scores2,pval<=0.05)
str(sig.sp.sc)
row.names(sig.sp.sc)=NULL
sig.sp.sc$Orders=as.character(sig.sp.sc$Orders)

#DEPRICATED Rename orders with codes for pretty plotting 
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Caulobacterales'] <- 'Ca'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Chroococcales'] <- 'Ch'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Corynebacteriales'] <- 'Co'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Gemmatales'] <- 'Ge'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Gemmatimonadales'] <- 'Gd'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Lactobacillales'] <- 'La'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Micrococcales'] <- 'Mi'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Micromonosporales'] <- 'Ms'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Myxococcales'] <- 'My'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Nostocales'] <- 'No'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Oscillatoriales'] <- 'Os'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Pirellulales'] <- 'Pi'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Planctomycetales'] <- 'Pl'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Pleurocapsales'] <- 'Pe'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Propionibacteriales'] <- 'Pr'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Pseudonocardiales'] <- 'Ps'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Sporadotrichida'] <- 'Sp'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Streptomycetales'] <- 'St'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Streptosporangiales'] <- 'Ss'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Synechococcales'] <- 'Sy'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Thalassiosirales'] <- 'Th'
sig.sp.sc$Orders[sig.sp.sc$Orders == 'Xanthomonadales'] <- 'Xa'

#Add phylum/class column (check every time this is rerun) for plotting by color
sig.sp.sc=sig.sp.sc%>%
  mutate(Phylum=c("Alphaproteobacteria",
                  "Cyanobacteria",
                  "Actinobacteria",
                  "Planctomycetota",
                  "Gemmatimonadota",
                  "Firmicutes",
                  "Actinobacteria",
                  "Actinobacteria",
                  "Deltaproteobacteria",
                  "Cyanobacteria",
                  "Cyanobacteria",
                  "Planctomycetota",
                  "Planctomycetota",
                  "Cyanobacteria",
                  "Actinobacteria",
                  "Actinobacteria",
                  "Actinobacteria",
                  "Actinobacteria",
                  "Cyanobacteria",
                  "Heterokonta",
                  "Gammaproteobacteria"))
#Make hull vectors and hull data frame for plotting polygons
age.1 = data.score[(data.score$grp)$Age == "1",][chull(data.score[(data.score$grp)$Age=="1",c("NMDS1","NMDS2")]),]
age.2=data.score[(data.score$grp)$Age == "3",][chull(data.score[(data.score$grp)$Age== "3",c("NMDS1","NMDS2")]),]
age.3=data.score[(data.score$grp)$Age == "5",][chull(data.score[(data.score$grp)$Age=="5",c("NMDS1","NMDS2")]),]
hull.data=rbind(age.1,age.2,age.3)
str(hull.data)

#Plot 
nmds3=ggplot() +
  geom_polygon(data=hull.data,
               aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
               alpha=0.4)+
  geom_polygon(data=data.score,
               aes(x=NMDS1,y=NMDS2, group=grp$Mat),
               size=1,
               linetype=3, 
               fill=NA,
               colour="#333333")+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  scale_color_manual(labels=c("09:00","21:00","09:00"),values=c("#90a295","#455765","#cd9fb2"),name="Time") +
  scale_fill_manual(values=c("1"="#90a295","3"="#455765","5"="#cd9fb2")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  #coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  annotate(geom="text",x=1.1,y=-1.5, label="Stress = 0.125", color="black",size=5)+
  annotate(geom="text",x=0,y=-1, label="3", color="#333333",size=5)+
  annotate(geom="text",x=0.78,y=-0.42, label="1", color="#333333",size=5)+
  annotate(geom="text",x=-0.8,y=0.2, label="4", color="#333333",size=5)+
  annotate(geom="text",x=0.3,y=1.25, label="2", color="#333333",size=5)+
  scale_x_continuous(name="NMDS1")+
  scale_y_continuous(name="NMDS2")
nmds3
dev.print(png,fil="nmds_top.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
#Plot with vector overlays
nmds4=ggplot() +
  #geom_polygon(data=hull.data,
  #             aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
  #             alpha=0.4)+
 # geom_polygon(data=data.score,
  #             aes(x=NMDS1,y=NMDS2, group=grp$Mat),
  #             size=1,
  #             linetype=3, 
  #             fill=NA,
  #             colour="#333333")+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat,colour=grp$Age),size=6) +
  ggrepel::geom_text_repel(data=sig.sp.sc,aes(x=NMDS1,y=NMDS2,label=Orders),
                           cex=5,
                           direction="both",
                           segment.size=0.2,
                           max.overlaps=50,
                           box.padding = 0.25,
                           force=40,
                           alpha=0.8,
                           segment.alpha=0.8,
                           min.segment.length = 0)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.sp.sc, size =1, alpha = 0.5, colour = "grey30")+
  scale_color_manual(labels=c("09:00","21:00","09:00"),values=c("#90a295","#455765","#cd9fb2"),name="Time") +
  scale_fill_manual(values=c("1"="#90a295","3"="#455765","5"="#cd9fb2")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  #coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  scale_y_continuous(name="NMDS2")
nmds4

#Color by phylum instead of having labels
nmds4=ggplot() +
  #geom_polygon(data=hull.data,
  #             aes(x=NMDS1,y=NMDS2,fill=grp$Age,group=grp$Age),
  #             alpha=0.4)+
  # geom_polygon(data=data.score,
  #             aes(x=NMDS1,y=NMDS2, group=grp$Mat),
  #             size=1,
  #             linetype=3, 
  #             fill=NA,
  #             colour="#333333")+
  geom_point(data=data.score,aes(x=NMDS1,y=NMDS2, shape=grp$Mat),size=6) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, color=Phylum), 
               data = sig.sp.sc, size =1)+
  scale_color_manual(values=c("#8acd65","#fae655","#47948b","#3d0b51","#3e6f8b","#dee14e","#3d6389","#56ac82","#432a73"), name="Phylum/Class") +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  #coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  scale_y_continuous(name="")+
  annotate(geom="text",x=1.1,y=-1.5, label="Stress = 0.125", color="black",size=5)
nmds4
dev.print(png,fil="nmds_bottom.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

library(patchwork)
dna_nmds_arrange=nmds3 + nmds4
dna_nmds_arrange=dna_nmds_arrange + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face='bold'),
        legend.position = "right")
dna_nmds_arrange=dna_nmds_arrange+plot_layout(guides="collect")
dna_nmds_arrange
dev.print(png,fil="nmds_stack.png",type="quartz",antialias="default",width=11,height=5.5,units="in",res=1300)


#Try using a jaccard index on raw data to look at turnover vs. nestedness
##Wrangle the data
turn_vs_nest_data=diel%>%
  select(-c(percent,taxon_id))%>%
  separate(taxon_name,c("Kingdom","Phylum","Class","Order"),";")%>%
  select(-c("Kingdom","Phylum","Class"))%>%
  filter(!is.na(Order))
str(turn_vs_nest_data)
turn_vs_nest_data$Order=as.factor(turn_vs_nest_data$Order)
turn_vs_nest_data2=turn_vs_nest_data%>%
  filter(Order != " Above Order" & Order != " <1%")
#Pivot
turn_vs_nest_wide=pivot_wider(turn_vs_nest_data2,names_from=Order,values_from=reads)
#Replace NA with 0
turn_vs_nest_wide[is.na(turn_vs_nest_wide)] <- 0
turn_vs_nest_wide$file=as.character(turn_vs_nest_wide$file)
#Add time and mat columns
turn_vs_nest_wide=turn_vs_nest_wide %>%
  mutate(Time=case_when(
    endsWith(file,"1") ~ "1",
    endsWith(file,"3") ~ "3",
    endsWith(file,"5") ~ "5",
  )) %>%
  mutate(Mat=case_when(
    str_detect(file,"D1") ~ "1",
    str_detect(file,"D2") ~ "2",
    str_detect(file,"D3") ~ "3",
    str_detect(file,"D4") ~ "4",))

#Make new dataframes for nmds
tt_dat=turn_vs_nest_wide[,2:75]
tt_dat[tt_dat > 0] <- 1
tt_grp=turn_vs_nest_wide[,76:77]
library(betapart)
#Run betapart with sorenson index
diel_betapart=beta.pair(tt_dat,index.family="sorensen")

#Extract turnover component
diel_betapart_turn=as.dist(diel_betapart$beta.sim)

#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(diel_betapart_turn, tt_grp$Time)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.2283

#Mat
dispr2 <- vegan::betadisper(diel_betapart_turn, tt_grp$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.6926
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(diel_betapart_turn~tt_grp$Time+tt_grp$Mat,permutations=999, method="sorenson", by="margin")
permanova #neither significant

#Extract nestedness component
diel_betapart_nest=as.dist(diel_betapart$beta.sne)

#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(diel_betapart_nest, tt_grp$Time)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.2283

#Mat
dispr2 <- vegan::betadisper(diel_betapart_nest, tt_grp$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.6926
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(diel_betapart_nest~tt_grp$Time+tt_grp$Mat,permutations=999, method="sorenson", by="margin")
permanova #Time significant

#Pairwise perm on signigificant Time
pair_perm=pairwiseAdonis::pairwise.adonis(diel_betapart_nest, factors=tt_grp$Time, sim.method="sorenson")
pair_perm #1v3

#Examine NMDS of this s
time_soren_nest_nmds=metaMDS(diel_betapart_nest,distance="sorenson",weakties=F)
time_soren_nest_nmds
stressplot(time_soren_nest_nmds,
           las=1,
           pch=19,
           cex=1,
           lwd=3)
#Fails, try PCA
time_soren_nest_pca=prcomp(diel_betapart_nest,scale=T)
#Scree
plot(time_soren_nest_pca) #Scree looks great 
library(ggfortify)
gg_time_soren=autoplot(time_soren_nest_pca,data=tt_grp,colour='Time',frame=T,shape='Mat',size=6)
gg_time_soren+theme_classic()+
  scale_color_manual(labels=c("09:00","21:00","09:00"),values=c("#90a295","#455765","#cd9fb2"),name="Time") +
  scale_fill_manual(values=c("1"="#90a295","3"="#455765","5"="#cd9fb2")) +
  scale_shape_manual(values=c(19,17,15,18),name="Mat")+
  #coord_equal()+
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  guides(fill=FALSE) +
  guides(shape=guide_legend(override.aes=list(size=5,color="#333333")))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  guides(size=F)
dev.print(png,fil="/Users/cissell/Desktop/sorenson_nest.png",type="quartz",antialias="default",width=6.5,
          height=5.5,units="in",res=1500)
#Total sorenson component
diel_betapart_sor=as.dist(diel_betapart$beta.sor)
#Extract mean and SD for pub
mean(diel_betapart_sor)#0.05
sd(diel_betapart_sor)#0.02
#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(diel_betapart_sor, tt_grp$Time)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.2283

#Mat
dispr2 <- vegan::betadisper(diel_betapart_sor, tt_grp$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.6926
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(diel_betapart_sor~tt_grp$Time+tt_grp$Mat,permutations=999, method="sorenson", by="margin")
permanova #Time significant

#Look for nestedness as evidence of succesioanl reset by time




#PCA Analysis
phylo_clr_rda<- ordinate(phylo_clr, method = "RDA", distance="euclidean")
phylo_clr_rda
#Scree plot from PCA
phyloseq::plot_scree(phylo_clr_rda) + 
  geom_bar(stat="identity", fill = "#455765") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")+
  theme_classic()+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1))
dev.print(png,fil="/Users/cissell/Desktop/scree.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

head(phylo_clr_rda$CA$eig)   

#PCA plot with circles around age (white background)
clr1 <- phylo_clr_rda$CA$eig[1] / sum(phylo_clr_rda$CA$eig)
clr2 <- phylo_clr_rda$CA$eig[2] / sum(phylo_clr_rda$CA$eig)
ord=plot_ordination(phylo_clr, phylo_clr_rda,
                     color = "Age", shape = "Mat") +
  theme_classic()
ord=ord+geom_point(size=5.9)+
  #stat_ellipse(aes(group = Age), linetype = 2, size=1.4, type="t")+
  coord_fixed(clr2 / clr1) +
  scale_color_manual(values=c("#90a295","#455765","#cd9fb2"))+
  #scale_color_manual(values=c("#acbad5","#e5d0b3"))+ #other two colors for RNA with 5 times
  #scale_color_manual(values=c("#87a8e2","#fea7ab","#fed1a2"))+ #Red, blue, yellow pastel
  scale_shape_manual(values=c(19,17,15,18))+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1))+
  guides(shape=guide_legend(override.aes=list(size=5)))
ord
dev.print(png,fil="ordination.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)


#Create distance matrix from clr
clr_dist_matrix <- phyloseq::distance(phylo_clr, method = "euclidean")
distance(phylo_clr, method = "euclidean")

#Check homogeneity of dispersion among groups
#Age
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(phylo_clr)$Age)
dispr
plot(dispr)
permutest(dispr,permutations=999) #p=0.2283

#Mat
dispr2 <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(phylo_clr)$Mat)
dispr2
plot(dispr2)
permutest(dispr2,permutations=999) #p=0.6926
#Assumption of homogeneity of dispersion holds true (non-sig p-value from permutest). Proceed with PERMANOVA

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(clr_dist_matrix ~sample_data(phylo_clr)$Age+sample_data(phylo_clr)$Mat,permutations=999, method="euclidean", by="margin")
permanova #Both significant: Age marginal = p 0.00085; Mat marginal = p 0.00570

clr_jaccard_matrix <- phyloseq::distance(phylo_clr, method = "jaccard")
#PERMANOVA on Age+Mat
permanova=vegan::adonis2(clr_jaccard_matrix ~sample_data(phylo_clr)$Age+sample_data(phylo_clr)$Mat,permutations=999, method="jaccard", by="margin")
permanova #Both significant: Age marginal = p 0.00085; Mat marginal = p 0.00570

#PERMANOVA on Age+Mat
permanova=vegan::adonis2(clr_dist_matrix ~sample_data(phylo_clr)$Age+sample_data(phylo_clr)$Mat,permutations=999, method="bray", by="margin")
permanova #Both significant: Age marginal = p 0.00085; Mat marginal = p 0.00570
#Beta-Diversity Network (look into hive-plots for other shit)

#Heatmap
heat_time=plot_heatmap(phylo_clr,sample.order=c("D1_1","D2_1","D3_1","D4_1","D1_3","D2_3","D3_3","D4_3","D1_5","D2_5","D3_5","D4_5"),low="#440154FF",high="#FDE725FF", taxa.label="Order")
heat_time=heat_time+theme_classic()
heat_time

heat_mat=plot_heatmap(phylo_clr,distance="euclidean",low="#440154FF",high="#FDE725FF", taxa.label="Order", method="PCoA")
heat_mat=heat_mat+theme_classic()
heat_mat

##Richness
#By mat
#We will add in Nonpareil Nd values, too.
nd=read.xls("/Users/cissell/Desktop/nd.xlsx")
est_rich=estimate_richness(phylo_object,measures=c("Chao1", "Shannon"))
est_rich$Nd=nd$Nd
est_rich$Age=nd$time
#Pivot!
est_rich_pivot=pivot_longer(data=est_rich, cols = c("Chao1","Shannon","Nd"), names_to="metric", values_to="value")
#Ggplot that s
rich_plot_all=ggplot(data=est_rich_pivot, aes(x=Age, y=value, group=Age)) +
  geom_boxplot(alpha=0.65,
               fill=c("#90a295","#455765","#cd9fb2",
                      "#90a295","#455765","#cd9fb2",
                      "#90a295","#455765","#cd9fb2"),
               outlier.shape=NA) +
  geom_point(aes(x=Age, y=value, colour=factor(Age)),size=4,position=position_dodge2(width=0.8,padding=0.1))+
  facet_wrap(.~metric, scales="free_y") +
  theme_classic() +
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt")) +
  theme(strip.text.x = element_text(size=15)) +
  theme(strip.background = element_blank()) +
  scale_color_manual(values=c("#90a295","#455765","#cd9fb2")) +
  theme(legend.position="none") +
  scale_x_continuous(breaks=c(1,3,5), name="Time Point",
                     labels=c("1"="09:00","2"="21:00","3"="09:00")) +
  scale_y_continuous(name="Diversity Metric Value")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  
rich_plot_all
dev.print(png,fil="/Users/cissell/Desktop/rich_plot.png",type="quartz",antialias="default",width=8.5,
          height=6.5,units="in",res=1300)
Fig_1=dna_nmds_arrange | rich_plot_all 
Fig_1=Fig_1+plot_layout(widths=c(1,1))
Fig_1 + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face='bold'))
dev.print(png,fil="/Users/cissell/Desktop/rich_plot.png",type="quartz",antialias="default",width=10,
          height=8.5,units="in",res=1500)
#Create matrix
rich=estimate_richness(phylo_object)
vec=c(1,2,3,1,2,3,1,2,3,1,2,3)
rich$Age = vec

#Test for differences
pairwise.wilcox.test(rich$Shannon,sample_data(phylo_object)$Age, p.adjust.method = "bonferroni") #no significant differences by age
pairwise.wilcox.test(rich$Shannon,sample_data(phylo_object)$Mat, p.adjust.method = "bonferroni") #no significant differenes by mat
pairwise.wilcox.test(rich$Chao1,sample_data(phylo_object)$Age, p.adjust.method = "bonferroni") #no significant differences by age
pairwise.wilcox.test(rich$Chao1,sample_data(phylo_object)$Mat, p.adjust.method = "bonferroni") #no significant differenes by mat
pairwise.wilcox.test(nd$Nd,nd$time,p.adjust.method = 'bonferroni') #no sif diff by age with nonpareil diversity index
#ggplot Heatmap
library(viridis)

heatmap(otu_table(phylo_clr),Rowv=NA,col=viridis(100), scale="column")

###Lets look at an association metric that replaces Pearson's correlation coefficient (proportianlity in propr package)
library(propr)

###DEPRICATED Use ALDEx2 to get differential abundance
library(ALDEx2)
conflict_prefer("mutate", "dplyr")
#set conditions list corresponding to order in otu_mat (input matrix)
conds=c("1","3","1","3","1","3","1","3")
#Create pairwise matrix, dropping all _5
otu_mat_2=subset(otu_mat,select=-c(D1_5,D2_5,D3_5,D4_5))
#Run ALDEx2
diff=aldex(otu_mat_2,conds,mc.samples=128,effect=TRUE,denom="all",verbose=TRUE,CI=T)
#Subset by significance (authors recommend effect size > |1| over p value)
sig_dif=subset(diff, effect >1 | effect < (-1))
#Create new column from rownames
sig_dif$row=row.names(sig_dif)
#Separate those into identifier columns
sig_dif=separate(sig_dif,row, c("Domain", "Phylum", "Class", "Order"),";")
#Sort by effect
sig_dif=sig_dif[order(sig_dif$effect),]
View(sig_dif)
#Reorder factor levels
sig_dif$Order=factor(sig_dif$Order, levels=(sig_dif$Order))

#plot sig dif
sig_dif_plot=ggplot(sig_dif,aes(x=effect,y=Order))+
  geom_point() +
  theme_classic()+
  xlab("Effect") +
  ylab("Order")+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.title.y = element_text(margin =
                                      margin(t = 0, r = 8, b = 0, l = 0)))+
  labs(y= "Order",
       fill= "Change")+
  scale_fill_manual(values=c("#90a295","#455765"))
sig_dif_plot
dev.print(png,fil="/Users/cissell/Desktop/diff_abund.png",type="quartz",antialias="default",width=6.5,
          height=8.5,units="in",res=1300)


#Run all again for 1 v 5 
conds2=c("1","5","1","5","1","5","1","5")
#Create pairwise matrix, dropping all _3
otu_mat_3=subset(otu_mat,select=-c(D1_3,D2_3,D3_3,D4_3))
#Run ALDEx2
diff2=aldex(otu_mat_3,conds2,mc.samples=128,effect=TRUE,denom="all",verbose=TRUE)
View(diff2) #None sig diff
#Subset by significance (authors recommend effect size > |1| over p value)
sig_dif2=subset(diff2, effect >1| effect < (-1))
#Create new column from rownames
sig_dif2$row=row.names(sig_dif2)
#Separate those into identifier columns
sig_dif2=separate(sig_dif2,row, c("Domain", "Phylum", "Class", "Order"),";")
#Sort by effect
sig_dif2=sig_dif2[order(sig_dif2$effect),]
View(sig_dif2)



#Run all again for 3 v 5
conds3=c("3","5","3","5","3","5","3","5")
#Create pairwise matrix, dropping all _1
otu_mat_4=subset(otu_mat,select=-c(D1_1,D2_1,D3_1,D4_1))
#Run ALDEx2
diff3=aldex(otu_mat_4,conds3,mc.samples=128,effect=TRUE,denom="all",verbose=TRUE)
View(sig_dif3)
#Subset by significance (authors recommend effect size > |1| over p value)
sig_dif3=subset(diff3, effect >1 | effect < (-1))
#Create new column from rownames
sig_dif3$row=row.names(sig_dif3)
#Separate those into identifier columns
sig_dif3=separate(sig_dif3,row, c("Domain", "Phylum", "Class", "Order"),";")
#Sort by effect
sig_dif3=sig_dif3[order(sig_dif3$effect),]
#Plot
sig_dif_plot2=ggplot(sig_dif3,aes(x=effect,y=reorder(Order,acsd(effect))))+
  geom_point() +
  theme_classic()+
  xlab("Effect") +
  ylab("Order")+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.title.y = element_text(margin =
                                      margin(t = 0, r = 8, b = 0, l = 0)))+
  labs(y= "Order",
       fill= "Change")+
  scale_fill_manual(values=c("#90a295","#455765"))
sig_dif_plot2

#We can also do a density plot of all of the pairwise distances (modified circa that nature diel coral pub that does it with dots instead, but mine is better)
#First get a dataframe with the distances
dist_data=as.data.frame(as.matrix(clr_dist_matrix))
#Write that out to a csv
write.csv(dist_data,file="/Users/cissell/Desktop/distance.csv")
#Now read that in
dist_xls=read.xls('/Users/cissell/Desktop/distance.xlsx')
#Now plot that in
#Density
library(ggpubr)
library(conflicted)
conflict_prefer("dplyr")
y_density=ggplot(data=dist_xls, aes(dist, fill=Comparison)) +
  geom_density(alpha=0.7,color="Black",size=1) +
  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2")) +
  scale_x_continuous(limits=c(1.5,4.7), name="")+
  ggpubr::rotate() +
  clean_theme() +
  theme(panel.background = element_blank())
y_density
#Scatter plot with points
scattered=ggplot()+
  geom_point(data=dist_xls, aes(x=x, y=dist, color=Comparison), size=5)+
  scale_color_manual(values=c("#90a295","#455765","#cd9fb2"))+
  theme_classic()+
  scale_y_continuous(limits=c(1.5,4.7), name="Euclidean dissimilarity")+
  scale_x_continuous(limits=c(1,66),name="", breaks=NULL)+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line.y = element_line(size=1.1),
        axis.line.x = element_line(size=0))
scattered
#Now arrange it into a masterpiece
ggarrange(scattered, y_density, 
          ncol = 2, nrow = 1,  align = "hv", 
          widths = c(2, 1),
          common.legend = TRUE)
dev.print(png,fil="distances.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
View(dist_xls)

#Anova for differeces 
aov1=aov(dist~Comparison,data=dist_xls)
summary(aov1)
#Now let's do that but with a violin plot, y'all
violin=ggplot(data=dist_xls, aes(x=Comparison, y=dist,fill=Comparison))+
  geom_violin(trim=FALSE,size=1,alpha=0.9)+
  geom_boxplot(width=0.1,fill="#333333",
               colour="#f2f2f2",
               outlier.shape=NA,
               size=0.9)+
  geom_jitter(shape=16, position=position_jitter(0.1),alpha=0.8)+
  theme_classic()+
  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2"))+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_y_continuous(name='Dissimilarity')+
  theme(legend.position="none")
violin
dev.print(png,fil="violin.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)


#PACKAGE CODYN for temporal community metrics
??codyn
?codyn::rate_change
?codyn::turnover
codyn_df=psmelt(phylo_clr)

#View(codyn_df)
codyn_df=codyn_df%>%
  select(-OTU,-Domain,-Phylum,-Class,-Sample)
str(codyn_df)
codyn_df$Age=as.numeric(codyn_df$Age)

#rate_change
rate_change_slopes=rate_change(df=codyn_df,time.var='Age',species.var='Order',abundance.var='Abundance',replicate.var='Mat')
rate_change_df=rate_change_interval(df=codyn_df,time.var='Age',species.var='Order',abundance.var='Abundance',replicate.var='Mat')

rate_plot=ggplot(rate_change_df, aes(interval, distance, color = Mat)) + 
  scale_color_manual(values=c("#add8e6","#cd9fb2","#455765","#90a295"))+
  geom_point(size=4,alpha=0.7) + 
  theme_classic() + 
  stat_smooth(method = "lm", se = F, size = 3)+
  scale_x_continuous(name="Interval",breaks=c(1,2))+
  scale_y_continuous(name="Euclidean Distance")+
  theme(axis.text = element_text(size=13,color="black"),
      axis.title=element_text(size=15,color="black"),
      axis.line = element_line(size=1.1))
rate_plot
dev.print(png,fil="rate.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
########################################
########################################
########################################
########################################
#####DIEL SEQUENCING DATA
#Now we will use TSS data for interpretale visualization but no statistics
##Read in data
diel=read.xls('/Users/cissell/Desktop/kaiju_summary.xlsx')

##Wrangle the data
otu=pivot_wider(diel,id_cols=taxon_name,names_from = file, values_from = percent)
#Replace NA with 0
otu[is.na(otu)] <- 0

#Remove Unclassified
otu=otu[-c(64),]

#Make OTU table for phyloseq
#Rownames have to be ID
as.data.frame(otu)
rownames(otu)<- otu$taxon_name


#Check counts
count <- colSums(otu[, c(2:ncol(otu))])
count

#Create taxonomy table
taxa<- otu %>% 
  select(taxon_name) %>% 
  separate(taxon_name, c("Domain", "Phylum", "Class", "Order"),
           ";")
str(taxa)



#Delete redundant 1st column
otu<- otu %>% 
  select(-taxon_name)
#Make factors not strings
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
str(taxa)
#Add first column of OTU table so OTU and Taxa tables match
taxa<- cbind(otu$taxon_name, taxa)
#Rename first column
colnames(taxa)[1]<- "taxon_name"
str(taxa)
#OTU IDs have to be rownames, just like OTU table
rownames(taxa)<- taxa$taxon_name

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
otu<- otu %>% 
  select(-taxon_name)

taxa<- taxa %>% 
  select(-taxon_name)

#Create metadata table
meta=read.xls('/Users/cissell/Desktop/metadata_table_diel.xlsx')
str(meta)
rownames(meta)<- meta$SampleID
meta<- meta %>% 
  select(-SampleID)
#Make proper levels
meta$Mat<- factor(meta$Mat, levels = c("D1", "D2", "D3","D4"))
meta$Age<- factor(meta$Age, levels = c("1", "3", "5"))

#Create matrices before creating phyloseq object
otu_mat<- as.matrix(otu)
tax_mat<- tax_table(as.matrix(taxa))

#Transform to phyloseq objects
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples<- sample_data(meta)

#And unite them into one object
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)
#Check and see if it makes sense
sample_sums(phylo_object)
sample_names(phylo_object)
rank_names(phylo_object) 
sample_variables(phylo_object) 
otu_table(phylo_object)[1:3, 1:2]
taxa_names(phylo_object)[1:5]
tax_table(phylo_object)
taxa_names(phylo_object)

#Heatmap
heat_time=plot_heatmap(phylo_object,sample.order=c("D1_1","D2_1","D3_1","D4_1","D1_3","D2_3","D3_3","D4_3","D1_5","D2_5","D3_5","D4_5"),low="#440154FF",high="#FDE725FF", taxa.label="Order")
heat_time=heat_time+theme_classic()
heat_time
dev.print(png,fil="heatmap.png",type="quartz",antialias="default",width=6.5,
          height=8.5,units="in",res=1300)

#Venn diagram for presence absence
library(MicEco)
#Time
venn_age=ps_venn(phylo_object,group='Age',type='percent',relative=FALSE,fraction=0,fill=c("#90a295","#455765","#cd9fb2"))
#Individual
venn_mat=ps_venn(phylo_object,group='Mat',type='percent',relative=FALSE,fraction=0,fill=c("#90a295","#455765","#cd9fb2","#fde992"))
#Combine
venns=ggarrange(venn_age,venn_mat,
                      ncol=2,nrow=1,align="hv",
                      labels=c("a","b"))
venns
annotate_figure(venns,
                left=text_grob("Phylum", color="Black", rot=90,size=15,vjust=(2.3),hjust=0.1))
dev.print(png,fil="venns.png",type="quartz",antialias="default",width=6.5,
          height=4.5,units="in",res=1300)

#Stacked Barplot
#First we need to backtransform to a dataframe to remove all the less abundant orders
barplot_df=psmelt(phylo_object)
str(barplot_df)
View(barplot_df)
#make order character, not factor
barplot_df$Order=as.character(barplot_df$Order)

#make new column with renamed low abundant taxa
conflict_prefer("mutate", "dplyr")
barplot_df=barplot_df%>%
  mutate(Order2 = replace(Order, Abundance < 1.7, " <1%" ))

#combine <1%
levels(barplot_df$Order2)[levels(barplot_df$Order2)=="<1%"] = "<1%"

#check unique names
unique(barplot_df$Order2)

#Reorder based on abundance
barplot_df$Order2=reorder(barplot_df$Order2,barplot_df$Abundance)

#Set order of samples and taxa for plot
level_order=factor(barplot_df$Sample, level=c("D1_1","D2_1","D3_1","D4_1","D1_3","D2_3","D3_3","D4_3","D1_5","D2_5","D3_5","D4_5"))

inverted_level_order=factor(barplot_df$Sample, level=c("D4_5","D3_5","D2_5","D1_5","D4_3","D3_3","D2_3","D1_3","D4_1","D3_1","D2_1","D1_1"))

taxa_order=factor(barplot_df$Order2, level=c(" <1%"," Above Order","Rhodospirillales","Chromatiales","Cytophagales","Desulfobacterales","Cellvibrionales","Rhizobiales","Alteromonadales","Nostocales","Flavobacteriales","Rhodobacterales","Oscillatoriales"))

#Now let's plot
stacked_bar=ggplot(barplot_df, aes(x=Abundance,y=inverted_level_order, fill=taxa_order)) +
  geom_bar(stat = "identity", width=1) +
  scale_fill_manual(values =#c("#e7ffac",
                              #"#bffcc6",
                             # "#97a2ff",
                             # "#aff8db",
                             # "#f6a6ff",
                             # "#ffffd1",
                             # "#6eb5ff",
                             # "#f3ffe3",
                             # "#ffabab",
                             # "#c4faf8",
                             # "#ff9cee",
                             # "#fbe4ff",
                             # "#b28dff")) +
                            c("#6eb5ff",
                               "#cceeff",
                               "#645e9d",
                               "#392b58",
                               "#2d0320",
                               "#d45087",
                               "#a05195",
                               "#f95d6a",
                               "#665191",
                               "#ff7c43",
                               "#6c969d",
                               "#ffa600",
                               "#003f5c")) +
  labs(y= "Relative abundance [%]",
       fill= "Order")+
  theme_classic()+
  xlab("Relative Abundance (%)") +
  ylab("Sample")+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.title.y = element_text(margin =
                                      margin(t = 0, r = 8, b = 0, l = 0)))
stacked_bar
#
dev.print(png,fil="stackedbar.png",type="quartz",antialias="default",width=5.5,
          height=5.5,units="in",res=1300)


#Now let's make a boxplot of viral abundance over time
#First backtransform clr dataset
viral_diel_df=psmelt(phylo_clr)
str(viral_diel_df)
viral_diel=ggplot(data=subset(viral_diel_df,Phylum==' Viruses'),
                  aes(x=Age, y=Abundance)) +
  geom_boxplot(aes(x=Age, y=Abundance))
viral_diel                  
             
#Let's make a dendrogram of the samples 
install.packages("ggdendro")
install.packages("dendextend")
install.packages("circlize")
library(dendextend)
library(ggdendro)
library(circlize)
hc=hclust(clr_dist_matrix, method = 'ward.D2')
dend=as.dendrogram(hc)
dend=dend %>%
  color_branches(k=3,col=c("#455765","#cd9fb2","#90a295")) %>%
  color_labels(k=3,col=c("#455765","#cd9fb2","#90a295")) %>%
  set("branches_lwd",c(6,6,6)) %>%
  set("labels_cex",1.6)
par(mar=rep(0,4))
circlize_dendrogram(dend,labels_track_height = NA, dend_track_height = 0.5,facing="outside")
dev.print(png,fil="dendrogram.png",type="quartz",antialias="default",width=5.5,
          height=5.5,units="in",res=1300)




#####PHLA LEVEL
#Now we will use TSS data for interpretale visualization but no statistics
##Read in data
diel=read.xls('/Users/cissell/Desktop/kaiju_summary_phyla.xlsx')

##Wrangle the data
otu=pivot_wider(diel,id_cols=taxon_name,names_from = file, values_from = percent)
#Replace NA with 0
otu[is.na(otu)] <- 0

#Remove Unclassified
otu=otu[-c(25),]

#Make OTU table for phyloseq
#Rownames have to be ID
as.data.frame(otu)
rownames(otu)<- otu$taxon_name


#Check counts
count <- colSums(otu[, c(2:ncol(otu))])
count

#Create taxonomy table
taxa<- otu %>% 
  select(taxon_name) %>% 
  separate(taxon_name, c("Domain", "Phylum"),
           ";")
str(taxa)



#Delete redundant 1st column
otu<- otu %>% 
  select(-taxon_name)
#Make factors not strings
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
str(taxa)
#Add first column of OTU table so OTU and Taxa tables match
taxa<- cbind(otu$taxon_name, taxa)
#Rename first column
colnames(taxa)[1]<- "taxon_name"
str(taxa)
#OTU IDs have to be rownames, just like OTU table
rownames(taxa)<- taxa$taxon_name

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
otu<- otu %>% 
  select(-taxon_name)

taxa<- taxa %>% 
  select(-taxon_name)

#Create metadata table
meta=read.xls('/Users/cissell/Desktop/metadata_table_diel.xlsx')
str(meta)
rownames(meta)<- meta$SampleID
meta<- meta %>% 
  select(-SampleID)
#Make proper levels
meta$Mat<- factor(meta$Mat, levels = c("D1", "D2", "D3","D4"))
meta$Age<- factor(meta$Age, levels = c("1", "3", "5"))

#Create matrices before creating phyloseq object
otu_mat<- as.matrix(otu)
tax_mat<- tax_table(as.matrix(taxa))

#Transform to phyloseq objects
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples<- sample_data(meta)

#And unite them into one object
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)
#Check and see if it makes sense
sample_sums(phylo_object)
sample_names(phylo_object)
rank_names(phylo_object) 
sample_variables(phylo_object) 
otu_table(phylo_object)[1:3, 1:2]
taxa_names(phylo_object)[1:5]
tax_table(phylo_object)
taxa_names(phylo_object)

#First we need to backtransform to a dataframe to remove all the less abundant orders
barplot_df=psmelt(phylo_object)
str(barplot_df)
View(barplot_df)
#make order character, not factor
barplot_df$Phylum=as.character(barplot_df$Phylum)
#make new column with renamed low abundant taxa
barplot_df=barplot_df%>%
  dplyr::mutate(Phylum2 = replace(Phylum, Abundance < 1, "<1%" ))

#combine <1%
levels(barplot_df$Phylum2)[levels(barplot_df$Phylum2)=="<1%"] = "<1%"

#check unique names
unique(barplot_df$Phylum2)
#Obtain mean % Abundances and SE within time
#View(barplot_df)
str(barplot_df)
conflict_prefer("summarise","dplyr")
phyla.sum <- barplot_df %>%
  group_by(Mat, Age, Phylum2) %>%
  dplyr::summarise(Abundsum = sum(Abundance))%>%
  group_by(Age, Phylum2) %>%
  dplyr::summarise(N = n(),
            mAbund = mean(Abundsum),
            seAbund = sd(Abundsum)/sqrt(N))
View(phyla.sum)
str(phyla.sum)
#Make barplot with phyla abundances and errors for abund >1%
#Inverse order
phyla.sum$Age=factor(phyla.sum$Age,levels=rev(levels(phyla.sum$Age)))
phyla.sum=phyla.sum %>% dplyr::arrange(phyla.sum$mAbund)
phyla.sum$Phylum2=factor(phyla.sum$Phylum2,levels=unique(phyla.sum$Phylum2))
View(phyla.sum)
#plot
mean_barplot_time=ggplot(data=subset(phyla.sum,mAbund>1 & Phylum2 !="Above Phyla"),aes(x=Phylum2,y=mAbund,fill=Age)) +
  geom_col(position="dodge",width=0.9) +
  geom_errorbar(aes(ymin=mAbund-seAbund, ymax=mAbund+seAbund),
                position=position_dodge(0.9),
                width=0,
                size=0.8)+
  coord_flip()+
  scale_fill_manual(labels=c("09:00","21:00","09:00"),values=c("#cd9fb2","#455765","#90a295"))+
  labs(x= "Relative abundance [%]",
       fill= "Time")+
  theme_classic()+
  ylab("Relative Abundance (%)") +
  xlab("")+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.title.y = element_text(margin =
                                      margin(t = 0, r = 8, b = 0, l = 0))) +
  guides(fill=guide_legend(reverse = TRUE))
mean_barplot_time

#Wilcox test for within time differnces in cyano and proteo
#1
target=c("Proteobacteria", "Cyanobacteria")
phyla.sum.1=dplyr::filter(barplot_df,Phylum2 %in% target & Age=="1")
str(phyla.sum.1)
stat.test1=phyla.sum.1%>%
  wilcox_test(Abundance~Phylum, paired=FALSE)
stat.test1 #p=0.0571 ; w=1
#3
phyla.sum.3=dplyr::filter(barplot_df,Phylum2 %in% target & Age=="3")
str(phyla.sum.3)
stat.test3=phyla.sum.3%>%
  wilcox_test(Abundance~Phylum, paired=FALSE)
stat.test3 #p=0.114 ; w=14
#5
phyla.sum.5=dplyr::filter(barplot_df,Phylum2 %in% target & Age=="5")
str(phyla.sum.5)
stat.test5=phyla.sum.5%>%
  wilcox_test(Abundance~Phylum, paired=FALSE)
stat.test5 #p=0.686 ; w=6

#all
phyla.sum.all=dplyr::filter(barplot_df,Phylum2 %in% target)
str(phyla.sum.all)
stat.testall=phyla.sum.all%>%
  wilcox_test(Abundance~Phylum, paired=FALSE)
stat.testall #p=0.686 ; w=6
#Obtain mean % Abundances and SE within mat
#View(barplot_df)
str(barplot_df)
conflict_prefer("summarise","dplyr")
phyla.sum.mat <- barplot_df %>%
  group_by(Mat, Age, Phylum2) %>%
  summarise(Abundsum = sum(Abundance))%>%
  group_by(Mat, Phylum2) %>%
  summarise(N = n(),
            mAbund = mean(Abundsum),
            seAbund = sd(Abundsum)/sqrt(N))
View(phyla.sum.mat)
str(phyla.sum.mat)
#Make barplot with phyla abundances and errors for abund >1%
#Inverse order
phyla.sum.mat$Mat=factor(phyla.sum.mat$Mat,levels=rev(levels(phyla.sum.mat$Mat)))
phyla.sum.mat=phyla.sum.mat %>% dplyr::arrange(phyla.sum.mat$mAbund)
phyla.sum.mat$Phylum2=factor(phyla.sum.mat$Phylum2,levels=unique(phyla.sum.mat$Phylum2))
View(phyla.sum.mat)
#plot
mean_barplot_mat=ggplot(data=subset(phyla.sum.mat,mAbund>1 & Phylum2 !="Above Phyla"),aes(x=Phylum2,y=mAbund,fill=Mat)) +
  geom_col(position="dodge",width=0.9) +
  geom_errorbar(aes(ymin=mAbund-seAbund, ymax=mAbund+seAbund),
                position=position_dodge(0.9),
                width=0,
                size=0.8)+
  coord_flip()+
  scale_fill_manual(values=c("#add8e6","#cd9fb2","#455765","#90a295"))+
  labs(x= "Relative abundance [%]",
       fill= "Mat")+
  theme_classic()+
  ylab("Relative Abundance (%)") +
  xlab("")+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.title.y = element_text(margin =
                                      margin(t = 0, r = 8, b = 0, l = 0))) +
  guides(fill=guide_legend(reverse = TRUE))
mean_barplot_mat 
dev.print(png,fil="mat_phyla.png",type="quartz",antialias="default",width=5.5,
          height=5.5,units="in",res=1300)
phyla_abund=ggarrange(mean_barplot_time,NULL,mean_barplot_mat,
          ncol=1,nrow=3,align="h",
          widths=c(1,-0.1,1),
          heights=c(1,-0.18,1.1),
          labels=c("a","","b"),hjust=1.5)
annotate_figure(phyla_abund,
                left=text_grob("Phylum", color="Black", rot=90,size=15,vjust=(2.3),hjust=0.1))
dev.print(png,fil="time_phyla.png",type="quartz",antialias="default",width=4.5,
          height=6.5,units="in",res=1300)

#Obtain overall means across all samples
phyla.sum.all <- barplot_df %>%
  group_by(Mat, Age, Phylum2) %>%
  summarise(Abundsum = sum(Abundance))%>%
  group_by(Phylum2) %>%
  summarise(N = n(),
            mAbund = mean(Abundsum),
            seAbund = sd(Abundsum)/sqrt(N))
View(phyla.sum.all)
str(phyla.sum)


###################
##################
#Cyano ANI Comparison
##################
###################
ANI=read.csv("/Users/cissell/Desktop/Cyano_ANI_all.csv")
mean(ANI$ANI)
sd(ANI$ANI)
#Plot
ANI_plot=ggplot(data=ANI, aes(x=Distance,y=ANI))+
  geom_point(size=4.5,color="#333333",alpha=0.9)+
  stat_summary(fun="mean",geom="line",size=2.5,color="#cd9fb2")+
  theme_classic()+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_x_discrete(breaks=c("1","2","3"),limits=c("1","2","3"),name="Relative Mat Linear Distance")
ANI_plot
dev.print(png,fil="ani.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)


###################
##################
#Cyano Growth Rate Comaprison
##################
###################
library(mgcv)
library(DHARMa)
library(GGally)
IREP=read.xls("/Users/cissell/Desktop/diel_irep_cyano.xlsx")
str(IREP)
IREP$Mat=as.factor(IREP$Mat)
#Run a GAM
gam=gam(IREP$irep~s(IREP$Time,k=3,bs="gp"),data=IREP, method="REML",
        family="gaussian",
        link="identiy")
summary(gam)
#Anova for main effects
anova(gam)
#Test for overdispersion
simulationOutput <- simulateResiduals(gam)
testDispersion(simulationOutput)
plot(simulationOutput)

#Plot
irep_plot=ggplot(data=IREP, aes(x=Time,y=irep, color=Mat))+
  stat_smooth(method="loess",formula=y~x,size=2, color="#333333",
              geom="smooth",
              fill="#455765")+
  geom_point(size=4.5,alpha=0.9, position = position_dodge(0.2))+
  geom_line(aes(color=Mat),linetype=2,alpha=0.8,size=1)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  xlab("Time Point")+
  ylab("iRep Value")+
  scale_color_manual(values=c("#9c8cdb","#cd9fb2","#455765","#90a295"))+
  scale_x_continuous(breaks=c(1,3,5))
irep_plot
dev.print(png,fil="irep.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
