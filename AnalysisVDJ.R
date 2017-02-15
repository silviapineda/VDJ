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
library(untb)
library(lme4)

setwd("/Users/Pinedasans/Documents/VDJ/")

load("/Users/Pinedasans/Data/VDJ/VDJ.Rdata")
load("/Users/Pinedasans/Data/VDJ/VDJ_downSampled.Rdata")
load("/Users/Pinedasans/Data/VDJ/VDJ_amplification_downSampled.Rdata")

##########################################################
####Calculate entropy with downsampling by individual ####
#########################################################

###Clonotype for gDNA 
data_clonesInference_gDNA<-data_clonesInference[which(data_clonesInference$amplification_template=="gDNA"),]
specimen_unique<-unique(data_clonesInference$specimen_label)

data_clone_specimen<-NULL
clone_entropy<-matrix(NA,length(specimen_unique),5)
simpson<-matrix(NA,length(specimen_unique),5)
for(j in 1:5){
  print(j)
  data_clone_dataset<-data_clonesInference_gDNA[which(data_clonesInference_gDNA$dataset==paste("ClonesInfered_downsampled_",j,".csv",sep="")),]
  for (i in 1:length(specimen_unique)){
    #print(i)
    data_specimen_unique<-data_clone_dataset[which(data_clone_dataset$specimen_label==specimen_unique[i]),]
    clone_entropy[i,j]<-entropy(as.numeric(data_specimen_unique[,"numberClone"]))
    simpson[i,j]<-simpson(table(as.numeric(data_specimen_unique[,"numberClone"])))
  }
}

clone_entropy_mean<-apply(clone_entropy,1,function(x) mean(x,na.rm=T))
clone_entropy_norm<-clone_entropy_mean/max(clone_entropy_mean,na.rm = T)
clonality<-(1-clone_entropy_norm)
names(clonality)<-specimen_unique

id.sample<-match(names(clonality),reads_clones_down_annot$specimen_id)
reads_clones_down_annot[id.sample,"clonality_gDNA"]<-clonality

###Clones longitudional
reads_clones_annot_down_Long<-reads_clones_down_annot[which(reads_clones_down_annot$clin!="AR" & reads_clones_down_annot$clin!="pre-AR"),]
clinLong<-factor(reads_clones_annot_down_Long$clin, levels=c("NP", "PNR", "PR"))

##gDNA
reads_clones_annot_Long_gDNA <- reads_clones_annot_down_Long[which(reads_clones_annot_down_Long$read_count_gDNA>=100),] ## 67 samples (Discard samples with less than 100 reads)
clinLong_gDNA <- clinLong[which(reads_clones_annot_down_Long$read_count_gDNA>=100)] #(Discard samples with less than 100 reads)

p1<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=clonality_gDNA, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality gDNA")

par(mfrow = c(2, 2))  
boxplot(reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)], 
        main="Time 0", ylab = "Clonality", col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,0.5))
fit = lm((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==0)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)], 
        main="Time 6",ylab = "Clonality", col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,0.5))
fit = lm((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)], ylab = "Clonality", 
        main="Time 24",col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,0.5))
fit = lm((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)])
anova(fit)

boxplot(reads_clones_annot_AR_gDNA$clonality_gDNA ~ clinAR_gDNA, main = "AR" , ylab = "Number of clones",col = c("goldenrod","firebrick3"),ylim=c(0,0.5))
fit = lm((reads_clones_annot_AR_gDNA$clonality_gDNA~ clinAR_gDNA))
anova(fit)




###Clonotype for cDNA 
data_clonesInference_cDNA<-data_clonesInference[which(data_clonesInference$amplification_template=="cDNA"),]
specimen_unique<-unique(data_clonesInference$specimen_label)

data_clone_specimen<-NULL
clone_entropy<-matrix(NA,length(specimen_unique),5)
simpson<-matrix(NA,length(specimen_unique),5)
for(j in 1:5){
  print(j)
  data_clone_dataset<-data_clonesInference_cDNA[which(data_clonesInference_cDNA$dataset==paste("ClonesInfered_downsampled_",j,".csv",sep="")),]
  for (i in 1:length(specimen_unique)){
    #print(i)
    data_specimen_unique<-data_clone_dataset[which(data_clone_dataset$specimen_label==specimen_unique[i]),]
    clone_entropy[i,j]<-entropy(as.numeric(data_specimen_unique[,"numberClone"]))
    simpson[i,j]<-simpson(table(as.numeric(data_specimen_unique[,"numberClone"])))
  }
}

clone_entropy_mean<-apply(clone_entropy,1,function(x) mean(x,na.rm=T))
clone_entropy_norm<-clone_entropy_mean/max(clone_entropy_mean,na.rm = T)
clonality<-(1-clone_entropy_norm)
names(clonality)<-specimen_unique

id.sample<-match(names(clonality),reads_clones_down_annot$specimen_id)
reads_clones_down_annot[id.sample,"clonality_cDNA"]<-clonality

###Clones longitudional
reads_clones_annot_down_Long<-reads_clones_down_annot[which(reads_clones_down_annot$clin!="AR" & reads_clones_down_annot$clin!="pre-AR"),]
clinLong<-factor(reads_clones_annot_down_Long$clin, levels=c("NP", "PNR", "PR"))

##cDNA
reads_clones_annot_Long_cDNA <- reads_clones_annot_down_Long[which(reads_clones_annot_down_Long$read_count_cDNA>=100),] ## 67 samples (Discard samples with less than 100 reads)
clinLong_cDNA <- clinLong[which(reads_clones_annot_down_Long$read_count_cDNA>=100)] #(Discard samples with less than 100 reads)

par(mfrow = c(1, 2))  
boxplot(reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==6)] ~ clinLong_gDNA[which(reads_clones_annot_Long_cDNA$time==6)], 
        main="Time 0", ylab = "Clonality", col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==6)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_cDNA$clonality_gDNA[which(reads_clones_annot_Long_cDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_cDNA$time==24)], ylab = "Clonality", 
        main="Time 24",col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,0.5))
fit = lm((reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==24)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==24)])
anova(fit)



##########################################################
####Calculate entropy with downsampling by sample ####
#########################################################

###Clonotype for gDNA 
data_clonesInference_gDNA<-data_clonesInference_amp[which(data_clonesInference_amp$amplification_template=="gDNA"),]
specimen_unique<-unique(data_clonesInference_amp$specimen_label)

data_clone_specimen<-NULL
clone_entropy<-matrix(NA,length(specimen_unique),5)
simpson<-matrix(NA,length(specimen_unique),5)
for(j in 1:5){
  print(j)
  data_clone_dataset<-data_clonesInference_gDNA[which(data_clonesInference_gDNA$dataset==paste("ClonesInfered_downsampled_byAmplification_",j,".csv",sep="")),]
  for (i in 1:length(specimen_unique)){
    #print(i)
    data_specimen_unique<-data_clone_dataset[which(data_clone_dataset$specimen_label==specimen_unique[i]),]
    clone_entropy[i,j]<-entropy(as.numeric(data_specimen_unique[,"numberClone"]))
    simpson[i,j]<-simpson(table(as.numeric(data_specimen_unique[,"numberClone"])))
  }
}

clone_entropy_mean<-apply(clone_entropy,1,function(x) mean(x,na.rm=T))
clone_entropy_norm<-clone_entropy_mean/max(clone_entropy_mean,na.rm = T)
clonality<-(1-clone_entropy_norm)
names(clonality)<-specimen_unique

id.sample<-match(names(clonality),reads_clones_down_annot_amp$specimen_id)
reads_clones_down_annot_amp[id.sample,"clonality_gDNA"]<-clonality

###Clones longitudional
reads_clones_annot_down_Long<-reads_clones_down_annot_amp[which(reads_clones_down_annot_amp$clin!="AR" & reads_clones_down_annot_amp$clin!="pre-AR"),]
clinLong<-factor(reads_clones_annot_down_Long$clin, levels=c("NP", "PNR", "PR"))

##gDNA
reads_clones_annot_Long_gDNA <- reads_clones_annot_down_Long[which(reads_clones_annot_down_Long$read_count_gDNA>=100),] ## 67 samples (Discard samples with less than 100 reads)
clinLong_gDNA <- clinLong[which(reads_clones_annot_down_Long$read_count_gDNA>=100)] #(Discard samples with less than 100 reads)

p1<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=clonality_gDNA, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality gDNA")

par(mfrow = c(2, 2))  
boxplot(reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)], 
        main="Time 0", ylab = "Clonality", col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==0)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)], 
        main="Time 6",ylab = "Clonality", col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)], ylab = "Clonality", 
        main="Time 24",col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_gDNA$clonality_gDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)])
anova(fit)

boxplot(reads_clones_annot_AR_gDNA$clonality_gDNA ~ clinAR_gDNA, main = "AR" , ylab = "Number of clones",col = c("goldenrod","firebrick3"),ylim=c(0,0.5))
fit = lm((reads_clones_annot_AR_gDNA$clonality_gDNA~ clinAR_gDNA))
anova(fit)

###Clonotype for cDNA 
data_clonesInference_cDNA<-data_clonesInference_amp[which(data_clonesInference_amp$amplification_template=="cDNA"),]
specimen_unique<-unique(data_clonesInference_amp$specimen_label)

data_clone_specimen<-NULL
clone_entropy<-matrix(NA,length(specimen_unique),5)
simpson<-matrix(NA,length(specimen_unique),5)
for(j in 1:5){
  print(j)
  data_clone_dataset<-data_clonesInference_cDNA[which(data_clonesInference_cDNA$dataset==paste("ClonesInfered_downsampled_byAmplification_",j,".csv",sep="")),]
  for (i in 1:length(specimen_unique)){
    #print(i)
    data_specimen_unique<-data_clone_dataset[which(data_clone_dataset$specimen_label==specimen_unique[i]),]
    clone_entropy[i,j]<-entropy(as.numeric(data_specimen_unique[,"numberClone"]))
    simpson[i,j]<-simpson(table(as.numeric(data_specimen_unique[,"numberClone"])))
  }
}

clone_entropy_mean<-apply(clone_entropy,1,function(x) mean(x,na.rm=T))
clone_entropy_norm<-clone_entropy_mean/max(clone_entropy_mean,na.rm = T)
clonality<-(1-clone_entropy_norm)
names(clonality)<-specimen_unique

id.sample<-match(names(clonality),reads_clones_down_annot_amp$specimen_id)
reads_clones_down_annot_amp[id.sample,"clonality_cDNA"]<-clonality

###Clones longitudional
reads_clones_annot_down_Long<-reads_clones_down_annot_amp[which(reads_clones_down_annot_amp$clin!="AR" & reads_clones_down_annot_amp$clin!="pre-AR"),]
clinLong<-factor(reads_clones_annot_down_Long$clin, levels=c("NP", "PNR", "PR"))

##cDNA
reads_clones_annot_Long_cDNA <- reads_clones_annot_down_Long[which(reads_clones_annot_down_Long$read_count_cDNA>=100),] ## 67 samples (Discard samples with less than 100 reads)
clinLong_cDNA <- clinLong[which(reads_clones_annot_down_Long$read_count_cDNA>=100)] #(Discard samples with less than 100 reads)

par(mfrow = c(1, 2))  
boxplot((reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==6)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==6)], 
        main="Time 6",ylab = "Clonality", col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==6)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==24)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==24)], ylab = "Clonality", 
        main="Time 24",col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==24)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==24)])
anova(fit)



############################################
####Calculate entropy without downsampling#
###########################################
##gDNA
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
specimen_unique<-unique(data_qc$specimen_label)

clone_entropy<-rep(NA,length(specimen_unique))
simpson<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
  clone_entropy[i]<-entropy(as.numeric(data_specimen_unique[,"igh_clone_id"]))
  simpson[i]<-simpson(table(as.numeric(data_specimen_unique[,"igh_clone_id"])))
}

clone_entropy_norm<-clone_entropy/max(clone_entropy,na.rm = T)
clonality<-(1-clone_entropy_norm)
names(clonality)<-specimen_unique


id.sample<-match(names(clonality),reads_clones_annot$specimen_id)
reads_clones_annot[id.sample,"clonality"]<-clonality

##Separate Long and AR
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]

clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))


####gDNA
reads_clones_annot_Long_gDNA<-reads_clones_annot_Long[which(reads_clones_annot_Long$gDNA_reads>=100),]
clinLong_gDNA<-factor(reads_clones_annot_Long_gDNA$clin, levels=c("NP", "PNR", "PR"))

p1<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=clonality, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality Long gDNA")


###Fitthing a longitudinal model
fm1 <- lmer(clonality~time + (time | Sample_id),data=reads_clones_annot_Long_gDNA)
fm2 <- lmer(clonality~ time*clin + (time | Sample_id),data=reads_clones_annot_Long_gDNA)
# also works with lme objects
anova(fm1, fm2) 

par(mfrow = c(2, 2))  
boxplot(reads_clones_annot_Long_gDNA$clonality[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main="Time 0")
boxplot(reads_clones_annot_Long_gDNA$clonality[which(reads_clones_annot_Long_gDNA$time==6)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main = "Time 6")
boxplot(reads_clones_annot_Long_gDNA$clonality[which(reads_clones_annot_Long_gDNA$time==24)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main="Time 24")

reads_clones_annot_AR_gDNA<-reads_clones_annot_AR[which(reads_clones_annot_AR$gDNA_reads>=100),]
clinAR_gDNA<-factor(reads_clones_annot_AR_gDNA$clin, levels=c("pre-AR","AR"))

p2<-ggplot(data=reads_clones_annot_AR_gDNA, aes(x=time, y=clonality, group=subject_id, shape=clinAR_gDNA, color=clinAR_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Clonality AR gDNA")

boxplot(reads_clones_annot_AR_gDNA$clonality ~ clinAR_gDNA , ylab = "Clonality", 
        col = c("goldenrod","firebrick3"))
fit = lm(reads_clones_annot_AR_gDNA$clonality ~ clinAR_gDNA)
anova(fit)

#####cDNA
reads_clones_annot_Long_cDNA<-reads_clones_annot_Long[which(reads_clones_annot_Long$cDNA_reads>=100),]
clinLong_cDNA<-factor(reads_clones_annot_Long_cDNA$clin, levels=c("NP", "PNR", "PR"))

p1<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=clonality, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() + 
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality Long cDNA")


reads_clones_annot_AR_cDNA<-reads_clones_annot_AR[which(reads_clones_annot_AR$cDNA_reads>=100),]
clinAR_cDNA<-factor(reads_clones_annot_AR_cDNA$clin, levels=c("pre-AR","AR"))

p2<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=clonality, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Clonality AR cDNA")




##cDNA
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
specimen_unique<-unique(data_qc$specimen_label)

clone_entropy<-rep(NA,length(specimen_unique))
simpson<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
  clone_entropy[i]<-entropy(as.numeric(data_specimen_unique[,"igh_clone_id"]))
  simpson[i]<-simpson(table(as.numeric(data_specimen_unique[,"igh_clone_id"])))
}

clone_entropy_norm<-clone_entropy/max(clone_entropy,na.rm = T)
clonality<-(1-clone_entropy_norm)
names(clonality)<-specimen_unique


id.sample<-match(names(clonality),reads_clones_annot$specimen_id)
reads_clones_annot[id.sample,"clonality_cDNA"]<-clonality

##Separate Long and AR
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]

clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))


####cDNA
reads_clones_annot_Long_cDNA<-reads_clones_annot_Long[which(reads_clones_annot_Long$cDNA_reads>=100),]
clinLong_cDNA<-factor(reads_clones_annot_Long_cDNA$clin, levels=c("NP", "PNR", "PR"))

par(mfrow = c(1, 2))  
boxplot(reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==6)] ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==6)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main = "Time 6")
boxplot(reads_clones_annot_Long_cDNA$clonality_cDNA[which(reads_clones_annot_Long_cDNA$time==24)] ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==24)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main="Time 24")

