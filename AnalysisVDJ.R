rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Analysis of the VDJ data
###         
###
### Author: Silvia Pineda
### Date: November, 2016
############################################################################################
library(entropy)
library(ggplot2)

setwd("/Users/Pinedasans/Documents/VDJ/")

load("/Users/Pinedasans/Data/VDJ/VDJ.Rdata")

samples_unique<-unique(data_qc$specimen_label)
clone_entropy<-NULL
for (i in 1:length(samples_unique)){
  data_sample_unique<-data_qc[which(data_qc$specimen_label==samples_unique[i]),"igh_clone_id"]
  clone_entropy[i]<-entropy(as.numeric(data_sample_unique))
}

clone_entropy_norm<-clone_entropy/max(clone_entropy)
clonality<-(1-clone_entropy_norm)

id.sample<-match(samples_unique,reads_clones_annot$specimen_id)
reads_clones_annot$clonality[id.sample]<-clonality

reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]

clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))

reads_clones_annot_Long_qc<- reads_clones_annot_Long[which(reads_clones_annot_Long$total_reads>=100),] ##  samples (Discard samples with less than 100 reads)
clinLong_qc <- clinLong[which(reads_clones_annot_Long$total_reads>=100)] #(Discard samples with less than 100 reads)

####Boxplot and ANOVA analysis for Longitudinal Data
###gDNA
reads_clones_annot_Long_gDNA<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$gDNA=="gDNA"),]
clinLong_gDNA<-clinLong_qc[which(reads_clones_annot_Long_qc$gDNA=="gDNA")]
boxplot(reads_clones_annot_Long_gDNA$clonality[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"))
boxplot(reads_clones_annot_Long_gDNA$clonality[which(reads_clones_annot_Long_gDNA$time==6)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"))
boxplot(reads_clones_annot_Long_gDNA$clonality[which(reads_clones_annot_Long_gDNA$time==24)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm(reads_clones_annot_Long_gDNA$clonality[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)])
anova(fit)

##AR 
reads_clones_annot_AR_qc<- reads_clones_annot_AR[which(reads_clones_annot_AR$total_reads>=100),] ##  samples (Discard samples with less than 100 reads)
clinAR_qc <- clinAR[which(reads_clones_annot_AR$total_reads>=100)] #(Discard samples with less than 100 reads)

boxplot(reads_clones_annot_AR_qc$clonality[which(reads_clones_annot_AR_qc$gDNA=="gDNA")] ~ clinAR_qc[which(reads_clones_annot_AR_qc$gDNA=="gDNA")] , ylab = "Clonality", col = c("goldenrod","firebrick3"))
fit = lm(reads_clones_annot_AR_qc$clonality[which(reads_clones_annot_AR_qc$gDNA=="gDNA")] ~ clinAR_qc[which(reads_clones_annot_AR_qc$gDNA=="gDNA")] )
anova(fit)

boxplot(reads_clones_annot_AR_qc$clonality[which(reads_clones_annot_AR_qc$cDNA=="cDNA")] ~ clinAR_qc[which(reads_clones_annot_AR_qc$cDNA=="cDNA")] , ylab = "Clonality", col = c("goldenrod","firebrick3"))
fit = lm(reads_clones_annot_AR_qc$clonality[which(reads_clones_annot_AR_qc$cDNA=="cDNA")] ~ clinAR_qc[which(reads_clones_annot_AR_qc$cDNA=="cDNA")] )
anova(fit)

###Longitudinal representation
ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=clonality, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality gDNA")


#AR
##gDNA
ggplot(data=reads_clones_annot_AR_qc[which(reads_clones_annot_AR_qc$gDNA=="gDNA"),] , aes(x=time, y=clonality, group=subject_id, shape=clinAR_qc[which(reads_clones_annot_AR_qc$gDNA=="gDNA")], color=clinAR_qc[which(reads_clones_annot_AR_qc$gDNA=="gDNA")])) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Clonality gDNA")

##cDNA
ggplot(data=reads_clones_annot_AR_qc[which(reads_clones_annot_AR_qc$cDNA=="cDNA"),] , aes(x=time, y=clonality, group=subject_id, shape=clinAR_qc[which(reads_clones_annot_AR_qc$cDNA=="cDNA")], color=clinAR_qc[which(reads_clones_annot_AR_qc$cDNA=="cDNA")])) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Clonality cDNA")