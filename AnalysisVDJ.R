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
library(sjPlot)

setwd("/Users/Pinedasans/VDJ/Results/")
load("/Users/Pinedasans/VDJ/Data/VDJ.Rdata")

####################################################
### Compare same individuals from gDNA and cDNA ###
###################################################


########################################
####Calculate repertoire diversity ####
#######################################

##########
## gDNA ##
##########

## 1. Clones by Cells
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]

##NP
specimen_unique_NP<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$clin=="NP")])
time<-unique(data_qc_gDNA$time[which(data_qc_gDNA$clin=="NP")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_NP)
  clones<-table(data_qc_gDNA$V_J_lenghCDR3_Clone_igh,data_qc_gDNA$specimen_label)
  
  for (i in 1:length(specimen_unique)){
    data_qc_gDNA_specimen<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_gDNA_specimen)[1]>=100){
      clones<-table(data_qc_gDNA_specimen$V_J_lenghCDR3_Clone_igh)
      write.table(clones,file=paste("NP-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[i],".txt",sep=""))
      tiff(paste("NP-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[i],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,1400),xlim=c(0,40),col=c("chartreuse4"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
      ggplot(data.matrix(table(clones)),aes())
      dev.off()
    }
  }
}
##PNR
specimen_unique_PNR<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$clin=="PNR")])
time<-unique(data_qc_gDNA$time[which(data_qc_gDNA$clin=="PNR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_PNR)
  for (i in 1:length(specimen_unique)){
    data_qc_gDNA_specimen<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_gDNA_specimen)[1]>=100){
      clones<-table(data_qc_gDNA_specimen$V_J_lenghCDR3_Clone_igh)
      data.matrix(table(clones))
      tiff(paste("PNR-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,1400),xlim=c(0,40),col=c("dodgerblue3"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}
##PR
specimen_unique_PR<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$clin=="PR")])
time<-unique(data_qc_gDNA$time[which(data_qc_gDNA$clin=="PR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_PR)
  for (i in 1:length(specimen_unique)){
    data_qc_gDNA_specimen<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_gDNA_specimen)[1]>=100){
      clones<-table(data_qc_gDNA_specimen$V_J_lenghCDR3_Clone_igh)
      tiff(paste("PR-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,1400),col=c("darkorange2"),pch=19,xlim=c(0,40),xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}
##AR
specimen_unique_AR<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$clin=="AR")])
time<-unique(data_qc_gDNA$time[which(data_qc_gDNA$clin=="AR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_AR)
  for (i in 1:length(specimen_unique)){
    data_qc_gDNA_specimen<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_gDNA_specimen)[1]>=100){
      clones<-table(data_qc_gDNA_specimen$V_J_lenghCDR3_Clone_igh)
      tiff(paste("AR-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,1400),col=c("darkorange2"),pch=19,xlim=c(0,40),xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}

##preAR
specimen_unique_preAR<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$clin=="pre-AR")])
time<-unique(data_qc_gDNA$time[which(data_qc_gDNA$clin=="pre-AR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_preAR)
  for (i in 1:length(specimen_unique)){
    data_qc_gDNA_specimen<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_gDNA_specimen)[1]>=100){
      clones<-table(data_qc_gDNA_specimen$V_J_lenghCDR3_Clone_igh)
      tiff(paste("pre-AR-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,1400),col=c("darkorange2"),pch=19,xlim=c(0,40),xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}


## 2. Diversity measures
#Species richness => Number of clones
#Entropy
#Clonality
#Simpson
specimen_unique<-unique(data_qc_gDNA$specimen_label)
entropy<-rep(NA,length(specimen_unique))
simpson<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  #print(i)
  data_specimen_unique<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
  clones_specimen<-data_specimen_unique[,"V_J_lenghCDR3_Clone_igh"]
  fi<-as.numeric(table(clones_specimen))/length(clones_specimen)
  hi<-fi*log2(fi)
  entropy[i]=-sum(hi)
  simpson[i]=sum(fi*fi)
  #entropy(table(clones_specimen)) returns the same result but by the fauls is the natural logarithm (log)
}


entropy_norm<-entropy/max(entropy,na.rm = T)
clonality<-(1-entropy_norm)
names(clonality)<-specimen_unique
diversity<-cbind(clonality,entropy,simpson)
write.csv(diversity,"diversity_gDNA.csv")


####Statistical Analysis
diversity<-read.table("/Users/Pinedasans/VDJ/Data/diversity.txt",sep="\t",header=T)

###Longitudinal Data
diversity_long<-diversity[which(diversity$clin=="NP" | diversity$clin=="PNR" | diversity$clin=="PR"),]
diversity_long_qc<-diversity_long[which(diversity_long$reads_gDNA>=100),]
clinLong<-factor(diversity_long_qc$clin, levels=c("NP", "PNR", "PR"))

p1<-ggplot(data=diversity_long_qc, aes(x=time, y=clones_gDNA, group=Sample_id, shape=clinLong, color=clinLong)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal gDNA")

diversity_long_qc$Sample_id_clin<-paste(diversity_long_qc$clin,"-",diversity_long_qc$Sample_id,sep="")
g1<-ggplot(diversity_long_qc[1:24,], aes(time, clones_gDNA)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_qc[25:48,], aes(time, clones_gDNA)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_qc[49:69,], aes(time, clones_gDNA)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")


###Fitthing a longitudinal model
fm1 <- lmer(diversity_long_qc$simpson_gDNA~time +clin+ (time | Sample_id),data=diversity_long_qc)
fm2 <- lmer(diversity_long_qc$simpson_gDNA~ time*clin + (time | Sample_id) ,data=diversity_long_qc)
ggplot(fortify(fm2), aes(time, diversity_long_qc$simpson_gDNA, color=clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
anova(fm1, fm2) 
sjp.glmer(fm2, type = "ri.pc",
          show.se = TRUE)
sjp.int(fm2)

###AR
diversity_AR<-diversity[which(diversity$clin=="AR" | diversity$clin=="pre-AR"),]
diversity_AR_qc<-diversity_AR[which(diversity_AR$reads_gDNA>=100),]
diversity_AR_qc$clin<-relevel(diversity_AR_qc$clin,ref="pre-AR")

fm1 <- lmer(diversity_AR_qc$entropy_gDNA~1+ (1 | Sample_id),data=diversity_AR_qc)
fm2 <- lmer(diversity_AR_qc$entropy_gDNA~ clin + (1 | Sample_id) ,data=diversity_AR_qc)
ggplot(fortify(fm2), aes(clin, diversity_AR_qc$entropy_gDNA, color=clin)) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
anova(fm1, fm2) 


##########
## cDNA ##
##########

## 1. Clones by Cells
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

## 2. Diversity measures
#Species richness => Number of clones
#Entropy
#Clonality
#Simpson
specimen_unique<-unique(data_qc_cDNA$specimen_label)
entropy<-rep(NA,length(specimen_unique))
simpson<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  #print(i)
  data_specimen_unique<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
  clones_specimen<-data_specimen_unique[,"V_J_lenghCDR3_Clone_igh"]
  fi<-as.numeric(table(clones_specimen))/length(clones_specimen)
  hi<-fi*log2(fi)
  entropy[i]=-sum(hi)
  simpson[i]=sum(fi*fi)
  #entropy(table(clones_specimen)) returns the same result but by the fauls is the natural logarithm (log)
}

entropy_norm<-entropy/max(entropy,na.rm = T)
clonality<-(1-entropy_norm)
names(clonality)<-specimen_unique
diversity<-cbind(clonality,entropy,simpson)
write.csv(diversity,"diversity_cDNA.csv")


####Statistical Analysis
diversity<-read.table("/Users/Pinedasans/VDJ/Data/diversity.txt",sep="\t",header=T)

###Longitudinal Data
diversity_long<-diversity[which(diversity$clin=="NP" | diversity$clin=="PNR" | diversity$clin=="PR"),]
diversity_long_qc<-diversity_long[which(diversity_long$reads_cDNA>=100),]
clinLong<-factor(diversity_long_qc$clin, levels=c("NP", "PNR", "PR"))

p1<-ggplot(data=diversity_long_qc, aes(x=time, y=clones_gDNA, group=Sample_id, shape=clinLong, color=clinLong)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal gDNA")

diversity_long_qc$Sample_id_clin<-paste(diversity_long_qc$clin,"-",diversity_long_qc$Sample_id,sep="")
ggplot(diversity_long_qc, aes(time, clones_cDNA)) + geom_point() + facet_grid(clin~Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")

###Fitthing a longitudinal model
fm1 <- lmer(diversity_long_qc$entropy_cDNA~time +clin+ (time | Sample_id),data=diversity_long_qc)
fm2 <- lmer(diversity_long_qc$entropy_cDNA~ time*clin + (time | Sample_id) ,data=diversity_long_qc)
ggplot(fortify(fm2), aes(time, diversity_long_qc$entropy_cDNA, color=clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
anova(fm1, fm2) 
sjp.int(fm2)

###AR
diversity_AR<-diversity[which(diversity$clin=="AR" | diversity$clin=="pre-AR"),]
diversity_AR_qc<-diversity_AR[which(diversity_AR$reads_cDNA>=100),]
diversity_AR_qc$clin<-relevel(diversity_AR_qc$clin,ref="pre-AR")

fm1 <- lmer(diversity_AR_qc$entropy_cDNA~1+ (1 | Sample_id),data=diversity_AR_qc)
fm2 <- lmer(diversity_AR_qc$entropy_cDNA~ clin + (1 | Sample_id) ,data=diversity_AR_qc)
ggplot(fortify(fm2), aes(clin, diversity_AR_qc$entropy_cDNA, color=clin)) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
anova(fm1, fm2) 
