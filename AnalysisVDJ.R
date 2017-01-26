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

samples_unique<-unique(data_qc$sample_id)
clone_entropy<-NULL
for (i in 1:length(samples_unique)){
  data_sample_unique<-data_qc[which(data_qc$sample_id==samples_unique[i]),"igh_clone_id"]
  clone_entropy[i]<-entropy(as.numeric(data_sample_unique))
}

clone_entropy_norm<-clone_entropy/max(clone_entropy)
clonality<-(1-clone_entropy_norm)

id.sample<-match(samples_unique,reads_clones_annot$Sample_id)
reads_clones_annot$clonality[id.sample]<-clonality

reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]

clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))

##Delete sample 28 for misclassification
reads_clones_annot_Long_2<-reads_clones_annot_Long[which(reads_clones_annot_Long$subject_id!=28),]

reads_clones_annot_Long_clones<- reads_clones_annot_Long_2[which(reads_clones_annot_Long_2$total_reads>=100),] ##  samples (Discard samples with less than 100 reads)
clinLong_clone <- clinLong[which(reads_clones_annot_Long_2$total_reads>=100)] #(Discard samples with less than 100 reads)


boxplot(reads_clones_annot_Long_clones$clonality ~ clinLong_clone, ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm(reads_clones_annot_Long_clones$clonality ~ clinLong_clone)
anova(fit)



###This will be only if we do per sample
reads_clones_annot_AR_clones<- reads_clones_annot_AR[which(reads_clones_annot_AR$total_reads>=100),] ##  samples (Discard samples with less than 100 reads)
clinAR_clone <- clinAR[which(reads_clones_annot_AR$total_reads>=100)] #(Discard samples with less than 100 reads)


boxplot(reads_clones_annot_AR_clones$total_reads ~ clinAR_clone, ylab = "Clonality", 
        col = c("goldenrod","firebrick3"))
fit = lm(reads_clones_annot_Long_clones$clonality ~ clinLong_clone)
anova(fit)




###Clones by individuals using all reads
ggplot(data=reads_clones_annot_Long_clones, aes(x=time, y=clonality, group=subject_id, shape=clinLong_clone, color=clinLong_clone)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality per Individual cDNA+gDNA")

boxplot((reads_clones_annot_Long_clones$clonality[which(reads_clones_annot_Long_clones$time==24)]) ~ clinLong_clone[which(reads_clones_annot_Long_clones$time==24)], ylab = "Clonality", 
            col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm(reads_clones_annot_Long_clones$clonality[which(reads_clones_annot_Long_clones$time==6)] ~ clinLong[which(reads_clones_annot_Long_clones$time==6)])
anova(fit)


