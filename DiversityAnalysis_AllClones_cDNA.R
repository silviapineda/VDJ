
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
library(gridExtra)
library(grid)
library(lattice)

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

###Number and mean by isotypes
table(table(data_cDNA_long$specimen_label)!=0) #56 samples
dim(data_cDNA_long)[1]/56
dim(data_cDNA_long[which(data_cDNA_long$isotype=="IGHE"),])[1]/56


###########################
## 2. Diversity measures###
###########################
############
## cDNA ### 
###########
specimen_unique<-unique(data_cDNA$specimen_label)
entropy_unmapped<-NULL
entropy_IGHA<-NULL
entropy_IGHD<-NULL
entropy_IGHE<-NULL
entropy_IGHG<-NULL
entropy_IGHM<-NULL
entropy_naive<-NULL
entropy_memory<-NULL
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_cDNA[which(data_cDNA$specimen_label==specimen_unique[i]),]
  clones_specimen_unmapped<-data_specimen_unique[which(data_specimen_unique$isotype==""),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHA<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHA"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHD<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHD"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHE<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHE"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHG<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHG"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHM<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHM"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_naive<-data_specimen_unique[which(data_specimen_unique$IGHM_naive_memory=="naive"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_memory<-data_specimen_unique[which(data_specimen_unique$IGHM_naive_memory=="memory"),"V_J_lenghCDR3_CloneId"]
  
  fi_unmapped<-as.numeric(table(clones_specimen_unmapped))/length(clones_specimen_unmapped)
  fi_IGHA<-as.numeric(table(clones_specimen_IGHA))/length(clones_specimen_IGHA)
  fi_IGHD<-as.numeric(table(clones_specimen_IGHD))/length(clones_specimen_IGHD)
  fi_IGHE<-as.numeric(table(clones_specimen_IGHE))/length(clones_specimen_IGHE)
  fi_IGHG<-as.numeric(table(clones_specimen_IGHG))/length(clones_specimen_IGHG)
  fi_IGHM<-as.numeric(table(clones_specimen_IGHM))/length(clones_specimen_IGHM)
  fi_naive<-as.numeric(table(clones_specimen_naive))/length(clones_specimen_naive)
  fi_memory<-as.numeric(table(clones_specimen_memory))/length(clones_specimen_memory)
  
  hi_unmapped<-fi_unmapped*log2(fi_unmapped)
  hi_IGHA<-fi_IGHA*log2(fi_IGHA)
  hi_IGHD<-fi_IGHD*log2(fi_IGHD)
  hi_IGHE<-fi_IGHE*log2(fi_IGHE)
  hi_IGHG<-fi_IGHG*log2(fi_IGHG)
  hi_IGHM<-fi_IGHM*log2(fi_IGHM)
  hi_naive<-fi_naive*log2(fi_naive)
  hi_memory<-fi_memory*log2(fi_memory)
  
  entropy_unmapped[i]=-sum(hi_unmapped)
  entropy_IGHA[i]=-sum(hi_IGHA)
  entropy_IGHD[i]=-sum(hi_IGHD)
  entropy_IGHE[i]=-sum(hi_IGHE)
  entropy_IGHG[i]=-sum(hi_IGHG)
  entropy_IGHM[i]=-sum(hi_IGHM)
  entropy_naive[i]=-sum(hi_naive)
  entropy_memory[i]=-sum(hi_memory)
  
}

#entropy_norm<-entropy/max(entropy,na.rm = T)
#clonality<-(1-entropy_norm)
#names(clonality)<-specimen_unique

diversity<-cbind(entropy_unmapped,entropy_IGHA,entropy_IGHD,entropy_IGHE,entropy_IGHG,entropy_IGHM,entropy_naive,entropy_memory)
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
diversity_long_cDNA$clin<-factor(diversity_long_cDNA$clin)
###Barplots
#reads
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
diversity_long_cDNA_IGHA<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHA>100),]
tiff("barplot_reads_IGHA_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHA$IGHA_isotypes[which(diversity_long_cDNA_IGHA$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==6)]],
        names.arg = diversity_long_cDNA_IGHA$Individual.id[which(diversity_long_cDNA_IGHA$time2==6)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHA)"),ylim = c(0,40000),cex.axis=0.8)
#legend(0.2, 40000, legend=levels(diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHA$IGHA_isotypes[which(diversity_long_cDNA_IGHA$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==24)]],
        names.arg = diversity_long_cDNA_IGHA$Individual.id[which(diversity_long_cDNA_IGHA$time2==24)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHA)"),ylim = c(0,40000),cex.axis=0.8)
dev.off()

diversity_long_cDNA_IGHD<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHD>100),]
tiff("barplot_reads_IGHD_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHD$IGHD_isotypes[which(diversity_long_cDNA_IGHD$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==6)]],
        names.arg = diversity_long_cDNA_IGHD$Individual.id[which(diversity_long_cDNA_IGHD$time2==6)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHD)"),ylim = c(0,70000),cex.axis=0.8)
barplot(diversity_long_cDNA_IGHD$IGHD_isotypes[which(diversity_long_cDNA_IGHD$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==24)]],
        names.arg = diversity_long_cDNA_IGHD$Individual.id[which(diversity_long_cDNA_IGHD$time2==24)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHD)"),ylim = c(0,70000),cex.axis=0.8)
#legend(0.2, 75000, legend=levels(diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==6)]),col=COLOR,pch=15, cex=1)
dev.off()

diversity_long_cDNA_IGHG<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHG>100),]
tiff("barplot_reads_IGHG_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHG$IGHG_isotypes[which(diversity_long_cDNA_IGHG$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==6)]],
        names.arg = diversity_long_cDNA_IGHG$Individual.id[which(diversity_long_cDNA_IGHG$time2==6)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHG)"),ylim = c(0,30000),cex.axis=0.8)
#legend(0.2, 30000, legend=levels(diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHG$IGHG_isotypes[which(diversity_long_cDNA_IGHG$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==24)]],
        names.arg = diversity_long_cDNA_IGHG$Individual.id[which(diversity_long_cDNA_IGHG$time2==24)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHG)"),ylim = c(0,30000),cex.axis=0.8)
dev.off()

diversity_long_cDNA_IGHM<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHM>100),]
tiff("barplot_reads_IGHM_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHM$IGHM_isotypes[which(diversity_long_cDNA_IGHM$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==6)]],
        names.arg = diversity_long_cDNA_IGHM$Individual.id[which(diversity_long_cDNA_IGHM$time2==6)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHM)"),ylim = c(0,60000),cex.axis=0.8)
#legend(0.2, 60000, legend=levels(diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHM$IGHM_isotypes[which(diversity_long_cDNA_IGHM$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==24)]],
        names.arg = diversity_long_cDNA_IGHM$Individual.id[which(diversity_long_cDNA_IGHM$time2==24)],
        cex.names=0.9,las=2,ylab = c("Reads (cDNA-IGHM)"),ylim = c(0,60000),cex.axis=0.8)
dev.off()

#clones
diversity_long_cDNA_IGHA<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHA>100),]
tiff("barplot_clones_IGHA_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHA$clones_IGHA[which(diversity_long_cDNA_IGHA$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==6)]],
        names.arg = diversity_long_cDNA_IGHA$Individual.id[which(diversity_long_cDNA_IGHA$time2==6)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHA)"),ylim = c(0,10000),cex.axis=0.8)
#legend(0.2, 40000, legend=levels(diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHA$clones_IGHA[which(diversity_long_cDNA_IGHA$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==24)]],
        names.arg = diversity_long_cDNA_IGHA$Individual.id[which(diversity_long_cDNA_IGHA$time2==24)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHA)"),ylim = c(0,10000),cex.axis=0.8)
dev.off()

diversity_long_cDNA_IGHD<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHD>100),]
tiff("barplot_clones_IGHD_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHD$clones_IGHD[which(diversity_long_cDNA_IGHD$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==6)]],
        names.arg = diversity_long_cDNA_IGHD$Individual.id[which(diversity_long_cDNA_IGHD$time2==6)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHD)"),ylim = c(0,50000),cex.axis=0.8)
barplot(diversity_long_cDNA_IGHD$clones_IGHD[which(diversity_long_cDNA_IGHD$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==24)]],
        names.arg = diversity_long_cDNA_IGHD$Individual.id[which(diversity_long_cDNA_IGHD$time2==24)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHD)"),ylim = c(0,50000),cex.axis=0.8)
#legend(0.2, 75000, legend=levels(diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==6)]),col=COLOR,pch=15, cex=1)
dev.off()

diversity_long_cDNA_IGHG<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHG>100),]
tiff("barplot_clones_IGHG_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHG$clones_IGHG[which(diversity_long_cDNA_IGHG$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==6)]],
        names.arg = diversity_long_cDNA_IGHG$Individual.id[which(diversity_long_cDNA_IGHG$time2==6)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHG)"),ylim = c(0,10000),cex.axis=0.7)
#legend(0.2, 30000, legend=levels(diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHG$clones_IGHG[which(diversity_long_cDNA_IGHG$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==24)]],
        names.arg = diversity_long_cDNA_IGHG$Individual.id[which(diversity_long_cDNA_IGHG$time2==24)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHG)"),ylim = c(0,10000),cex.axis=0.7)
dev.off()

diversity_long_cDNA_IGHM<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHM>100),]
tiff("barplot_clones_IGHM_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHM$clones_IGHM[which(diversity_long_cDNA_IGHM$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==6)]],
        names.arg = diversity_long_cDNA_IGHM$Individual.id[which(diversity_long_cDNA_IGHM$time2==6)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHM)"),ylim = c(0,40000),cex.axis=0.7)
#legend(0.2, 60000, legend=levels(diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHM$clones_IGHM[which(diversity_long_cDNA_IGHM$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==24)]],
        names.arg = diversity_long_cDNA_IGHM$Individual.id[which(diversity_long_cDNA_IGHM$time2==24)],
        cex.names=0.9,las=2,ylab = c("Clones (cDNA-IGHM)"),ylim = c(0,40000),cex.axis=0.7)
dev.off()

#SHM
diversity_long_cDNA_IGHA<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHA>100),]
diversity_long_cDNA_IGHA$SHM_IGHA_byclones<-diversity_long_cDNA_IGHA$SHM_IGHA/diversity_long_cDNA_IGHA$clones_IGHA

tiff("barplot_SHM_IGHA_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHA$SHM_IGHA_byclones[which(diversity_long_cDNA_IGHA$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==6)]],
        names.arg = diversity_long_cDNA_IGHA$Individual.id[which(diversity_long_cDNA_IGHA$time2==6)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHA)"),ylim = c(0,5))
#legend(0.2, 5, legend=levels(diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHA$SHM_IGHA_byclones[which(diversity_long_cDNA_IGHA$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHA$clin[which(diversity_long_cDNA_IGHA$time2==24)]],
        names.arg = diversity_long_cDNA_IGHA$Individual.id[which(diversity_long_cDNA_IGHA$time2==24)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHA)"),ylim = c(0,5))
dev.off()

diversity_long_cDNA_IGHD<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHD>100),]
diversity_long_cDNA_IGHD$SHM_IGHD_byclones<-diversity_long_cDNA_IGHD$SHM_IGHD/diversity_long_cDNA_IGHD$clones_IGHD
tiff("barplot_SHM_IGHD_cDNA.tiff",res=300,w=2500,h=3000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHD$SHM_IGHD_byclones[which(diversity_long_cDNA_IGHD$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==6)]],
        names.arg = diversity_long_cDNA_IGHD$Individual.id[which(diversity_long_cDNA_IGHD$time2==6)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHD)"),ylim = c(0,0.04))
barplot(diversity_long_cDNA_IGHD$SHM_IGHD_byclones[which(diversity_long_cDNA_IGHD$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==24)]],
        names.arg = diversity_long_cDNA_IGHD$Individual.id[which(diversity_long_cDNA_IGHD$time2==24)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHD)"),ylim = c(0,0.04))
#legend(0.2, 75000, legend=levels(diversity_long_cDNA_IGHD$clin[which(diversity_long_cDNA_IGHD$time2==6)]),col=COLOR,pch=15, cex=1)
dev.off()

diversity_long_cDNA_IGHG<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHG>100),]
diversity_long_cDNA_IGHG$SHM_IGHG_byclones<-diversity_long_cDNA_IGHG$SHM_IGHG/diversity_long_cDNA_IGHG$clones_IGHG
tiff("barplot_SHM_IGHG_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHG$SHM_IGHG_byclones[which(diversity_long_cDNA_IGHG$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==6)]],
        names.arg = diversity_long_cDNA_IGHG$Individual.id[which(diversity_long_cDNA_IGHG$time2==6)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHG)"),ylim = c(0,7))
#legend(0.2, 30000, legend=levels(diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHG$SHM_IGHG_byclones[which(diversity_long_cDNA_IGHG$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHG$clin[which(diversity_long_cDNA_IGHG$time2==24)]],
        names.arg = diversity_long_cDNA_IGHG$Individual.id[which(diversity_long_cDNA_IGHG$time2==24)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHG)"),ylim = c(0,7))
dev.off()

diversity_long_cDNA_IGHM<-diversity_long_cDNA[which(diversity_long_cDNA$clones_IGHM>100),]
diversity_long_cDNA_IGHM$SHM_IGHM_byclones<-diversity_long_cDNA_IGHM$SHM_IGHM/diversity_long_cDNA_IGHM$clones_IGHM
tiff("barplot_SHMs_IGHM_cDNA.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(diversity_long_cDNA_IGHM$SHM_IGHM_byclones[which(diversity_long_cDNA_IGHM$time2==6)],
        col = COLOR[diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==6)]],
        names.arg = diversity_long_cDNA_IGHM$Individual.id[which(diversity_long_cDNA_IGHM$time2==6)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHM)"),ylim = c(0,0.04))
#legend(0.2, 60000, legend=levels(diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA_IGHM$SHM_IGHM_byclones[which(diversity_long_cDNA_IGHM$time2==24)],
        col = COLOR[diversity_long_cDNA_IGHM$clin[which(diversity_long_cDNA_IGHM$time2==24)]],
        names.arg = diversity_long_cDNA_IGHM$Individual.id[which(diversity_long_cDNA_IGHM$time2==24)],
        cex.names=0.9,las=2,ylab = c("SHM (cDNA-IGHM)"),ylim = c(0,0.04))
dev.off()
###proportion of isotypes by NP-PNR-RR



####Filter out all those subjects with clones per isotype <=100
##Analysis by clones
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

j=1
PNR_6_clones<-NULL
PNR_24_clones<-NULL
PR_6_clones<-NULL
PR_24_clones<-NULL
for (i in c("clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")){
  print(i)
   assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(i,colnames(diversity_long_cDNA))]>100),])
  tiff(paste0("boxplot_clones_cDNA_",i,".tiff"),h=1800,w=3000,res=300)
  p2 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),], 
            aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)]), 
                get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) + 
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clinical outcome", y = "Number of clones") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none") + ylim(0,30000) +  theme(text = element_text(size=15)) 
  
   p3 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),], 
            aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)]), 
                get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) +
     geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clinical outcome", y = "Number of clones") + 
     stat_summary(fun.data=data_summary)  + theme(legend.position="none") + ylim(0,30000) + theme(text = element_text(size=15)) 
   
   grid.arrange(p2,p3,ncol=3)
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
names(PR_6_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")
names(PR_24_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")
names(PNR_6_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")
names(PNR_24_clones)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")

##Entropy
j=1
PNR_6_entropy<-NULL
PNR_24_entropy<-NULL
PR_6_entropy<-NULL
PR_24_entropy<-NULL
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")
for (i in c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM","entropy_naive","entropy_memory")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(clones[j],colnames(diversity_long_cDNA))]>100),])
  tiff(paste0("boxplot_entropy_cDNA_",i,".tiff"),h=2000,w=1800,res=300) 
  p2 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) + 
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clinical outcome", y = "Shannon entropy") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none")  + theme(text = element_text(size=15)) + ylim(7,15)
  
  p3 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) +
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clinical outcome", y = "Shannon entropy") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none")  + theme(text = element_text(size=15)) + ylim(7,15)
  
  grid.arrange(p2,p2,p3,ncol=3)
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
names(PR_6_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM","entropy_naive","entropy_memory")
names(PR_24_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM","entropy_naive","entropy_memory")
names(PNR_6_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM","entropy_naive","entropy_memory")
names(PNR_24_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM","entropy_naive","entropy_memory")

##SHM
j=1
PNR_6_SHM<-NULL
PNR_24_SHM<-NULL
PR_6_SHM<-NULL
PR_24_SHM<-NULL
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")
for (i in c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM","SHM_naive","SHM_memory")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(clones[j],colnames(diversity_long_cDNA))]>100),])
  tiff(paste0("boxplot_SHM_cDNA_",i,".tiff"),h=2000,w=1800,res=300)
  p2 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==6)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==6),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) + 
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clinical outcome", y = "SHM") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none")  + theme(text = element_text(size=15)) 
  
  p3 = ggplot(get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),], 
              aes(factor(get(paste0("diversity_long_cDNA_",i))$clin[which(get(paste0("diversity_long_cDNA_",i))$time2==24)]), 
                  get(paste0("diversity_long_cDNA_",i))[which(get(paste0("diversity_long_cDNA_",i))$time2==24),match(i,colnames(get(paste0("diversity_long_cDNA_",i))))],fill=clin)) +
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clinical outcome", y = "SHM") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none")  + theme(text = element_text(size=15)) 
  
  grid.arrange(p2,p3,ncol=2)
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
names(PR_6_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM","SHM_naive","SHM_memory")
names(PR_24_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM","SHM_naive","SHM_memory")
names(PNR_6_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM","SHM_naive","SHM_memory")
names(PNR_24_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM","SHM_naive","SHM_memory")

PR_6<-cbind(PR_6_clones,PR_6_entropy,PR_6_SHM)
PR_24<-cbind(PR_24_clones,PR_24_entropy,PR_24_SHM)
PNR_6<-cbind(PNR_6_clones,PNR_6_entropy,PNR_6_SHM)
PNR_24<-cbind(PNR_24_clones,PNR_24_entropy,PNR_24_SHM)
result<-cbind(PR_6,PR_24,PNR_6,PNR_24)
rownames(result)<-c("UNMAPPED","IGHA","IGHD","IGHG","IGHM")
write.csv(result,file="ResultsBoxplots_cDNA.csv")


########################################
###  Fitting a longitudinal model  ####
#######################################

##Clones
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")
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
  
  #tiff(paste0("plot_lmer_clones_cDNA_",i,".tiff"),h=1200,w=1400,res=300)
  p <- ggplot(fm_full, aes(x = time, y = get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1], colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=2) +
  geom_smooth(method="lm",size=1) +
  labs(x = "Time (months)",y = "Number of clones") + theme_bw() + theme_light()
  
  
  print(p)
  #dev.off()
}
names(p_value_clones)<-clones

##Entropy
clones<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM","clones_naive","clones_memory")
j<-1
p_value_entropy<-NULL
for (i in c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM","entropy_naive","entropy_memory")){
  assign(paste0("diversity_long_cDNA_",i),diversity_long_cDNA[which(diversity_long_cDNA[match(clones[j],colnames(diversity_long_cDNA))]>100),])
  fm_null <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin + time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  fm_full <- lmer(get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1] ~ 
                    clin*time + (1 | Sample_id),data=get(paste0("diversity_long_cDNA_",i)),REML = F)
  
  p_value_entropy[j]<-anova(fm_full, fm_null)[2,8] 
  j=j+1
  
  tiff(paste0("plot_lmer_entropy_cDNA_",i,".tiff"),h=1200,w=1400,res=300)
  p <- ggplot(fm_full, aes(x = time, y = get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1], colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=1.2) +
    geom_smooth(method="lm",size=0.8) +
    labs(x = "Time (months)",y = "Shannon entropy") + theme_bw() + theme_light()
  
  print(p)
  dev.off()
}
names(p_value_entropy)<-c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM","entropy_naive","entropy_memory")

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
  
  
  #tiff(paste0("plot_lmer_SHM_cDNA_",i,".tiff"),h=1200,w=1400,res=300)
  p <- ggplot(fm_full, aes(x = time, y = get(paste0("diversity_long_cDNA_",i))[match(i,colnames(get(paste0("diversity_long_cDNA_",i))))][,1], colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=3) +
    geom_smooth(method="lm",size=1.5) +
    labs(x = "time (months)",y = "SHM") 
  
  print(p)
  #dev.off()
}
names(p_value_SHM)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")


###Proportion of Naive and Memory by clin
diversity_long_cDNA$prop<-diversity_long_cDNA$clones_naive/diversity_long_cDNA$clones_IGHM

boxplot(diversity_long_cDNA$prop~factor(diversity_long_cDNA$clin))
summary(lm(diversity_long_cDNA$prop~factor(diversity_long_cDNA$clin)))

##Proportion of each subtype
diversity_long_cDNA$prop_IGHA<-diversity_long_cDNA$clones_IGHA/diversity_long_cDNA$clones_cDNA
boxplot(diversity_long_cDNA$prop_IGHA~factor(diversity_long_cDNA$clin))
summary(lm(diversity_long_cDNA$prop_IGHA~factor(diversity_long_cDNA$clin)))

diversity_long_cDNA$prop_IGHE<-diversity_long_cDNA$clones_IGHE/diversity_long_cDNA$clones_cDNA
boxplot(diversity_long_cDNA$prop_IGHE~factor(diversity_long_cDNA$clin))
summary(lm(diversity_long_cDNA$prop_IGHE~factor(diversity_long_cDNA$clin)))

diversity_long_cDNA$prop_IGHG<-diversity_long_cDNA$clones_IGHG/diversity_long_cDNA$clones_cDNA
boxplot(diversity_long_cDNA$prop_IGHG~factor(diversity_long_cDNA$clin))
summary(lm(diversity_long_cDNA$prop_IGHG~factor(diversity_long_cDNA$clin)))

diversity_long_cDNA$prop_IGHD<-diversity_long_cDNA$clones_IGHD/diversity_long_cDNA$clones_cDNA
boxplot(diversity_long_cDNA$prop_IGHD~factor(diversity_long_cDNA$clin))
summary(lm(diversity_long_cDNA$prop_IGHD~factor(diversity_long_cDNA$clin)))

diversity_long_cDNA$prop_IGHM<-diversity_long_cDNA$clones_IGHM/diversity_long_cDNA$clones_cDNA
boxplot(diversity_long_cDNA$prop_IGHM~factor(diversity_long_cDNA$clin))
summary(lm(diversity_long_cDNA$prop_IGHM~factor(diversity_long_cDNA$clin)))
