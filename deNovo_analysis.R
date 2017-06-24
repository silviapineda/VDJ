rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies. deNovo analysis
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Analysis of the VDJ data
###         
###
### Author: Silvia Pineda
### Date: June, 2017
############################################################################################
library(entropy)
library(ggplot2)
library(untb)
library(lme4)
library(sjPlot)
library("caroline")
library("effects")
library("RColorBrewer")

setwd("/Users/Pinedasans/VDJ/Results/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")

##############################
### gDNA longitudinal data ###
##############################

data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
#data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

### longitudinal
data_qc_gDNA_long<-data_qc_gDNA[which(data_qc_gDNA$clin!="AR" & data_qc_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_gDNA>100),]

id<-match(data_qc_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_gDNA_long_qc<-data_qc_gDNA_long[which(is.na(id)==F),]
data_qc_gDNA_long_qc$specimen_label<-factor(data_qc_gDNA_long_qc$specimen_label)

###Delete the samples that do not have three time points

###7_s63
data_qc_gDNA_long_qc2<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label!="7_S63"),]
##38021
data_qc_gDNA_long_qc3<-data_qc_gDNA_long_qc2[which(data_qc_gDNA_long_qc2$sample_id!="38021" & data_qc_gDNA_long_qc2$sample_id!="12008"
                                                   & data_qc_gDNA_long_qc2$sample_id!="38020" & data_qc_gDNA_long_qc2$sample_id!="41006"
                                                   & data_qc_gDNA_long_qc2$sample_id!="89002" & data_qc_gDNA_long_qc2$sample_id!="12001"
                                                   & data_qc_gDNA_long_qc2$sample_id!="89015" & data_qc_gDNA_long_qc2$sample_id!="301003"
                                                   & data_qc_gDNA_long_qc2$sample_id!="303001"),]

data_qc_gDNA_long_clean<-data_qc_gDNA_long_qc3

##Obtain deNovo number of clones
sample_id<-unique(data_qc_gDNA_long_clean$sample_id)
deNovo<-NULL
for (i in 1:length(sample_id)){
  print(i)
  clone_type_matrix<-as.data.frame.matrix(table(data_qc_gDNA_long_clean$V_J_lenghCDR3_Clone_igh[which(data_qc_gDNA_long_clean$sample_id==sample_id[i])],
                                      factor(data_qc_gDNA_long_clean$specimen_label[which(data_qc_gDNA_long_clean$sample_id==sample_id[i])])))
  clone_type_matrix$deNovo<-ifelse(clone_type_matrix[,1] == 0 & (clone_type_matrix[,2] !=0 | clone_type_matrix[,3] !=0),1,0)
  deNovo[i]<-sum(clone_type_matrix$deNovo)
}
names(deNovo)<-sample_id

###Read the clean annotated data
annot<-read.csv("/Users/Pinedasans/VDJ/Data/Annot_clean_long_gDNA.csv")
data_annot<-annot[match(names(deNovo),annot$Sample_id),c(1,2,4,6:24)]

#clin
summary(glm(deNovo~data_annot$clin))
boxplot(deNovo~data_annot$clin,col=c("chartreuse4", "dodgerblue3","darkorange2"))

##Multivariable analysis
summary(glm(deNovo~clin+immunosuppression+recipient.Race+GenderMismatch+
              Recipient.Age.when.had.Tx+Donor.Source+hla_mismatch,data=data_annot))

#immunosuppression
COLOR=brewer.pal(5,"Set1")
summary(glm(deNovo~data_annot$immunosuppression))
boxplot(deNovo~data_annot$immunosuppression,col=COLOR)

#Recipient Race
summary(glm(deNovo~data_annot$recipient.Race))
boxplot(deNovo~data_annot$recipient.Race,col=COLOR)

#Gender mismatch
summary(glm(deNovo~data_annot$GenderMismatch))
boxplot(deNovo~data_annot$GenderMismatch,col=COLOR)

#Age Recipient
summary(glm(deNovo~data_annot$Recipient.Age.when.had.Tx))
plot(deNovo~data_annot$Recipient.Age.when.had.Tx)

#Donor.Source
summary(glm(deNovo~data_annot$Donor.Source))
boxplot(deNovo~data_annot$Donor.Source,col=COLOR)

#hla_mismatch
summary(glm(deNovo~data_annot$hla_mismatch))
boxplot(deNovo~data_annot$hla_mismatch,col=COLOR)
