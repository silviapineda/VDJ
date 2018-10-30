rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies from Scott
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: September, 2018
############################################################################################
library(entropy)
library(ggplot2)
library(untb)
library(lme4)
library(caroline)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(lattice)
library(lmerTest)
library(finalfit)
library(dplyr)

COLOR=c("chartreuse4","dodgerblue3","darkorange2")


setwd("/Users/Pinedasans/VDJ/")
load("Data/diversity_long_gDNA.Rdata")

# Table 1 - Patient demographics by variable of interest ----
annot_qc_unique<-diversity_long_gDNA[which(duplicated(diversity_long_gDNA$Sample_id)==F),]
table(annot_qc_unique$clin)
explanatory = c("Recipient.Age.when.had.Tx", "immunosuppression","Donor.Source", "Recipient.Gender")
dependent = 'clin'
table<-annot_qc_unique %>%
  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)
tiff("Results/Table_demographic.tiff",res=300,w=4000,h=1200)
grid.table(table)
dev.off()

# Table 2 - Time by variable of interest ----
explanatory = c("time")
dependent = 'clin'
table<-diversity_long_gDNA %>%
  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)
tiff("Results/Table_time.tiff",res=300,w=4000,h=1200)
grid.table(table)
dev.off()



############################################
### Association with all the variables ####
###########################################
annot_qc_time0<-diversity_long_gDNA[which(diversity_long_gDNA$time==0),]
annot_qc_time6<-diversity_long_gDNA[which(diversity_long_gDNA$time>=6),]
#####Recipient Age
tiff("Results/boxplot_clones_RecAge.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(Recipient.Age.when.had.Tx,clones_gDNA,color=clin)) + 
  geom_point() + geom_smooth(method='lm') + labs(title="time <=0",x="Rec.Age", y = "Clones") +
  scale_color_manual(values=c("chartreuse4","dodgerblue3","darkorange2"))

p2 = ggplot(annot_qc_time6,aes(Recipient.Age.when.had.Tx,clones_gDNA,color=clin)) + 
  geom_point() + geom_smooth(method='lm') + labs(title="time >=6",x="Rec.Age", y = "Clones") +
  scale_color_manual(values=c("chartreuse4","dodgerblue3", "darkorange2"))
grid.arrange(p1,p2,ncol=2)
dev.off()

explanatory = c("clin","Recipient.Age.when.had.Tx")
dependent = 'clones_gDNA'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecAge_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecAge_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

#####Immunosupression
tiff("Results/boxplot_clones_Immuno.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(factor(immunosuppression),clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time <=0",x="Immunosuppression", y = "Clones")
p2 = ggplot(annot_qc_time6,aes(factor(immunosuppression),clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time >=6",x="Immunosuppression", y = "Clones")
grid.arrange(p1,p2,ncol=2)
dev.off()

explanatory = c("clin","immunosuppression")
dependent = 'clones_gDNA'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_Immuno_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_Immuno_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()


#####DonorType
tiff("Results/boxplot_clones_DonorType.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(factor(Donor.Source),clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time <=0",x="DonorType", y = "Clones")
p2 = ggplot(annot_qc_time6,aes(factor(Donor.Source),clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time >=6",x="DonorType", y = "Clones")
grid.arrange(p1,p2,ncol=2)
dev.off()

explanatory = c("clin","Donor.Source")
dependent = 'clones_gDNA'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_DonorType_0.tiff",res=300,w=4000,h=700)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_DonorType_6.tiff",res=300,w=4000,h=700)
grid.table(example_table)
dev.off()

#####RecipientSex
tiff("Results/boxplot_clones_RecipientSex.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(factor(Recipient.Gender),clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time <=0",x="RecipientSex", y = "Clones")
p2 = ggplot(annot_qc_time6,aes(factor(Recipient.Gender),clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time >=6",x="RecipientSex", y = "Clones")
grid.arrange(p1,p2,ncol=2)
dev.off()

explanatory = c("clin","Recipient.Gender")
dependent = 'clones_gDNA'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecipientSex_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecipientSex_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()


##########Considering all variables
explanatory = c("clin","Recipient.Age.when.had.Tx", "immunosuppression","Donor.Source", "Recipient.Gender")
dependent = 'clones_gDNA'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_All_0.tiff",res=300,w=4000,h=1000)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_All_6.tiff",res=300,w=4000,h=1000)
grid.table(example_table)
dev.off()



