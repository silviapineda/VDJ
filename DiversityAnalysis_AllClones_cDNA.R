
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
library("caroline")
library("effects")
library("RColorBrewer")

setwd("/Users/Pinedasans/VDJ/ResultsAllClones/")
load("/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")



########################################
####Calculate repertoire diversity ####
#######################################

##########
## cDNA ##
##########
data_cDNA<-data_merge[which(data_merge$amplification_template=="cDNA"),]
data_cDNA_long<-data_cDNA[which(data_cDNA$clin=="NP" | data_cDNA$clin=="PNR" | data_cDNA$clin=="PR"),]

###########################
## 2. Diversity measures###
###########################
############
## gDNA ### 
###########
specimen_unique<-unique(data_cDNA$specimen_label)
entropy_unmapped<-NULL
entropy_IGHA<-NULL
entropy_IGHD<-NULL
entropy_IGHE<-NULL
entropy_IGHG<-NULL
entropy_IGHM<-NULL
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_cDNA[which(data_cDNA$specimen_label==specimen_unique[i]),]
  clones_specimen_unmapped<-data_specimen_unique[which(data_specimen_unique$isotype==""),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHA<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHA"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHD<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHD"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHE<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHE"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHG<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHG"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHM<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHM"),"V_J_lenghCDR3_CloneId"]
  fi_unmapped<-as.numeric(table(clones_specimen_unmapped))/length(clones_specimen_unmapped)
  fi_IGHA<-as.numeric(table(clones_specimen_IGHA))/length(clones_specimen_IGHA)
  fi_IGHD<-as.numeric(table(clones_specimen_IGHD))/length(clones_specimen_IGHD)
  fi_IGHE<-as.numeric(table(clones_specimen_IGHE))/length(clones_specimen_IGHE)
  fi_IGHG<-as.numeric(table(clones_specimen_IGHG))/length(clones_specimen_IGHG)
  fi_IGHM<-as.numeric(table(clones_specimen_IGHM))/length(clones_specimen_IGHM)
  hi_unmapped<-fi_unmapped*log2(fi_unmapped)
  hi_IGHA<-fi_IGHA*log2(fi_IGHA)
  hi_IGHD<-fi_IGHD*log2(fi_IGHD)
  hi_IGHE<-fi_IGHE*log2(fi_IGHE)
  hi_IGHG<-fi_IGHG*log2(fi_IGHG)
  hi_IGHM<-fi_IGHM*log2(fi_IGHM)
  
  entropy_unmapped[i]=-sum(hi_unmapped)
  entropy_IGHA[i]=-sum(hi_IGHA)
  entropy_IGHD[i]=-sum(hi_IGHD)
  entropy_IGHE[i]=-sum(hi_IGHE)
  entropy_IGHG[i]=-sum(hi_IGHG)
  entropy_IGHM[i]=-sum(hi_IGHM)
}
#entropy_norm<-entropy/max(entropy,na.rm = T)
#clonality<-(1-entropy_norm)
#names(clonality)<-specimen_unique

diversity<-cbind(entropy,entropy_unmapped,entropy_IGHA,entropy_IGHD,entropy_IGHE,entropy_IGHG,entropy_IGHM)
rownames(diversity)<-specimen_unique
write.csv(diversity,"/Users/Pinedasans/VDJ/Data/diversity_AllClones_cDNA.csv")


########################################
####  Statistical Analysis  ############
#######################################

diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity_AllClones_cDNA.csv",header=T)

reads_clones_annot_cDNA<-reads_clones_annot[which(reads_clones_annot$cDNA=="cDNA"),]
id<-match(reads_clones_annot_cDNA$specimen_id,diversity$X)
diversity_reads_clones<-cbind(reads_clones_annot_cDNA,diversity[id,])

#############
### cDNA ####
#############

diversity_cDNA<-diversity_reads_clones[which(diversity_reads_clones$clones_cDNA>=1000),]

###Longitudinal Data
diversity_long_cDNA<-diversity_cDNA[which(diversity_cDNA$clin=="NP" | diversity_cDNA$clin=="PNR" | diversity_cDNA$clin=="PR"),]
diversity_long_cDNA$clin<-factor(diversity_long_cDNA$clin, levels=c("NP", "PNR", "PR"))
table(diversity_long_cDNA$clin[which(is.na(diversity_long_cDNA$reads_cDNA)==F)],diversity_long_cDNA$time[which(is.na(diversity_long_cDNA$reads_cDNA)==F)])

ggplot(data=diversity_long_cDNA, aes(x=time, y=clones_cDNA, group=Sample_id, shape=clin, color=clin)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal cDNA")


##################################
#####Analysis by time and clin ##
##################################

diversity_long_cDNA$time2<-replace(diversity_long_cDNA$time,diversity_long_cDNA$time==12,6)

###Barplots
#reads
for (i in c("UNMAPPED_isotypes","IGHA_isotypes","IGHD_isotypes","IGHE_isotypes","IGHG_isotypes","IGHM_isotypes")){
  tiff(paste0("barplot_reads_",i,"_cDNA.tiff"),res=300,w=3000,h=2500)
  par(mfrow=c(2,1))
  #barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),match(i,colnames(diversity_long_cDNA))],
   #     col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]],
    #    names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==0)],
     #   cex.names=0.8,las=2,ylab = c("Reads"))
  #legend(0.2, 250000, legend=levels(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]),col=COLOR,pch=15, cex=1)
  barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==6)],
        cex.names=0.8,las=2,ylab = c("Reads"))
  barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),match(i,colnames(diversity_long_cDNA))],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==24)],
        cex.names=0.8,las=2,ylab = c("Reads"))
  dev.off()
}

#clones
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
for (i in c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM")){
  tiff(paste0("barplot_",i,"_cDNA.tiff"),res=300,w=3000,h=2500)
  par(mfrow=c(2,1))
  #barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),match(i,colnames(diversity_long_cDNA))],
   #       col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]],
    #      names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==0)],
     #     cex.names=0.8,las=2,ylab = c("Reads"))
  #legend(0.2, 250000, legend=levels(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]),col=COLOR,pch=15, cex=1)
  barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))],
          col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]],
          names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==6)],
          cex.names=0.8,las=2,ylab = c("Clones"))
  barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),match(i,colnames(diversity_long_cDNA))],
          col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]],
          names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==24)],
          cex.names=0.8,las=2,ylab = c("Clones"))
  dev.off()
}

#SHM
diversity_long_cDNA$SHM_cDNA_byClones<-diversity_long_cDNA$SHM_cDNA/diversity_long_cDNA$clones_cDNA
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
for (i in c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHE","SHM_IGHG","SHM_IGHM")){
  tiff(paste0("barplot_",i,"_cDNA.tiff"),res=300,w=3000,h=2500)
  par(mfrow=c(2,1))
  #barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),match(i,colnames(diversity_long_cDNA))],
  #       col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]],
  #      names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==0)],
  #     cex.names=0.8,las=2,ylab = c("Reads"))
  #legend(0.2, 250000, legend=levels(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]),col=COLOR,pch=15, cex=1)
  barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))],
          col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]],
          names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==6)],
          cex.names=0.8,las=2,ylab = c("SHM"))
  barplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),match(i,colnames(diversity_long_cDNA))],
          col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]],
          names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==24)],
          cex.names=0.8,las=2,ylab = c("SHM"))
  dev.off()
}

####Filter out all those subjects with clones per isotype <=100
##Analysis by clones
j=1
PNR_6_clones<-NULL
PNR_24_clones<-NULL
PR_6_clones<-NULL
PR_24_clones<-NULL
for (i in c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(i,colnames(diversity_long_cDNA))]>100),])
  tiff(paste0("boxplot_clones_cDNA_",i,".tiff"),h=2000,w=1800,res=300)
  p2 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),], 
            aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)]), 
                get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
  p3 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),], 
            aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)]), 
                get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")

  multiplot(p2,p3)
  dev.off()
  PNR_6_clones[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))] ~ 
                  get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)])))[2,4]
  PR_6_clones[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))] ~ 
                   get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)])))[3,4]

  PNR_24_clones[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))] ~ 
                  get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)])))[2,4]
  PR_24_clones[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))] ~ 
                   get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)])))[3,4]

  j=j+1
}
names(PR_6_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
names(PR_24_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
names(PNR_6_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
names(PNR_24_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")

##Entropy
j=1
PNR_6_entropy<-NULL
PNR_24_entropy<-NULL
PR_6_entropy<-NULL
PR_24_entropy<-NULL
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
for (i in c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(clones[j],colnames(diversity_long_cDNA))]>100),])
  tiff(paste0("boxplot_entropy_cDNA_",i,".tiff"),h=2000,w=1800,res=300)
  p2 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) + 
    geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "entropy")
  p3 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) +
    geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "entropy")
  
  multiplot(p2,p3)
  dev.off()
  PNR_6_entropy[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))] ~ 
                               get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)])))[2,4]
  PR_6_entropy[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))] ~ 
                              get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)])))[3,4]
  
  PNR_24_entropy[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))] ~ 
                                get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)])))[2,4]
  PR_24_entropy[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))] ~ 
                               get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)])))[3,4]
  
  j=j+1
}
names(PR_6_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")
names(PR_24_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")
names(PNR_6_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")
names(PNR_24_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")

##SHM
j=1
PNR_6_SHM<-NULL
PNR_24_SHM<-NULL
PR_6_SHM<-NULL
PR_24_SHM<-NULL
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
for (i in c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(clones[j],colnames(diversity_long_cDNA))]>100),])
  tiff(paste0("boxplot_SHM_cDNA_",i,".tiff"),h=2000,w=1800,res=300)
  p2 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) + 
    geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "SHM")
  p3 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) +
    geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "SHM")
  
  multiplot(p2,p3)
  dev.off()
  PNR_6_SHM[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))] ~ 
                               get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)])))[2,4]
  PR_6_SHM[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(diversity_long_cDNA$time2==6),match(i,colnames(diversity_long_cDNA))] ~ 
                              get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)])))[3,4]
  
  PNR_24_SHM[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))] ~ 
                                get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)])))[2,4]
  PR_24_SHM[j]<-coef(summary(glm(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))] ~ 
                               get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)])))[3,4]
  
  j=j+1
}
names(PR_6_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")
names(PR_24_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")
names(PNR_6_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")
names(PNR_24_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")

PR_6<-cbind(PR_6_clones,PR_6_entropy,PR_6_SHM)
PR_24<-cbind(PR_24_clones,PR_24_entropy,PR_24_SHM)
PNR_6<-cbind(PNR_6_clones,PNR_6_entropy,PNR_6_SHM)
PNR_24<-cbind(PNR_24_clones,PNR_24_entropy,PNR_24_SHM)
result<-cbind(PR_6,PR_24,PNR_6,PNR_24)
rownames(result)<-c("UNMAPPED","IGHA","IGHD","IGHG","IGHM")
write.csv(result,file="ResultsBoxplots_cDNA.csv")


########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
j<-1
p_value_clones<-NULL
for (i in clones){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(i,colnames(diversity_long_cDNA))]>100),])
  fm_null <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin + time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  fm_full <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin*time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  
  p_value_clones[j]<-anova(fm_full, fm_null)[2,8] 
  j=j+1
  
  tiff(paste0("plot_lmer_clones_cDNA_",i,".tiff"),h=1700,w=2000,res=300)
  p <- ggplot(fm_full, aes(x = time, y = get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1], colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=3) +
    geom_smooth(method="lm",size=1.5) +
    labs(x = "time (months)",y = paste0("Clones (cDNA-",substr(i,9,13),")")) 
  
  print(p)
  dev.off()
}
names(p_value_clones)<-clones

##Entropy
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
j<-1
p_value_entropy<-NULL
for (i in c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(clones[j],colnames(diversity_long_cDNA))]>100),])
  fm_null <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin + time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  fm_full <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin*time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  
  p_value_entropy[j]<-anova(fm_full, fm_null)[2,8] 
  j=j+1
  
  tiff(paste0("plot_lmer_entropy_cDNA_",i,".tiff"),h=1700,w=2000,res=300)
  p <- ggplot(fm_full, aes(x = time, y = get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1], colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=3) +
    geom_smooth(method="lm",size=1.5) +
    labs(x = "time (months)",y = paste0("Entropy (cDNA-",substr(i,9,13),")")) 
  
  print(p)
  dev.off()
}
names(p_value_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")

##SHM
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
j<-1
p_value_SHM<-NULL
for (i in c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(clones[j],colnames(diversity_long_cDNA))]>100),])
  fm_null <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin + time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  fm_full <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin*time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  
  p_value_SHM[j]<-anova(fm_full, fm_null)[2,8] 
  j=j+1
  
  tiff(paste0("plot_lmer_SHM_cDNA_",i,".tiff"),h=1700,w=2000,res=300)
  p <- ggplot(fm_full, aes(x = time, y = get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1], colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=3) +
    geom_smooth(method="lm",size=1.5) +
    labs(x = "time (months)",y = "SHM") 
  
  print(p)
  dev.off()
}
names(p_value_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")







