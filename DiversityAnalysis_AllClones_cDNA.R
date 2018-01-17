
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
id<-match(reads_clones_annot_cDNA$specimen_id,rownames(diversity))
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
#clones
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
tiff("barplot_clones_time_cDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==0)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==0)],
        cex.names=0.8,las=2,ylim = c(0,70000),ylab = c("Clones"))
#legend(0, 80000, legend=levels(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==6)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==6)],
        cex.names=0.8,las=2,ylim = c(0,70000),ylab = c("Clones"))
barplot(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==24)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==24)],
        cex.names=0.8,las=2,ylim = c(0,70000),ylab = c("Clones"))
dev.off()

#reads
tiff("barplot_reads_time_cDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_cDNA$IGHD_isotypes[which(diversity_long_cDNA$time2==0)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==0)],
        cex.names=0.8,las=2,ylab = c("Reads"))
#legend(0.2, 250000, legend=levels(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA$IGHD_isotypes[which(diversity_long_cDNA$time2==6)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==6)],
        cex.names=0.8,las=2,ylab = c("Reads"))
barplot(diversity_long_cDNA$IGHD_isotypes[which(diversity_long_cDNA$time2==24)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==24)],
        cex.names=0.8,las=2,ylab = c("Reads"))
dev.off()


#SHM
diversity_long_cDNA$SHM_cDNA_byClones<-diversity_long_cDNA$SHM_cDNA/diversity_long_cDNA$clones_cDNA
tiff("barplot_SHM_byClones_cDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_cDNA$SHM_cDNA_byClones[which(diversity_long_cDNA$time2==0)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==0)],
        cex.names=0.8,las=2,ylim = c(0,0.2),ylab = c("SHM"))
legend(0, .2, legend=levels(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA$SHM_cDNA_byClones[which(diversity_long_cDNA$time2==6)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==6)],
        cex.names=0.8,las=2,ylim = c(0,0.2),ylab = c("SHM"))
barplot(diversity_long_cDNA$SHM_cDNA_byClones[which(diversity_long_cDNA$time2==24)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==24)],
        cex.names=0.8,las=2,ylim = c(0,0.2),ylab = c("SHM"))
dev.off()

###Sample 8 is an outlier at 6 and 24
diversity_long_gDNA_sample8<-diversity_long_gDNA[which(diversity_long_gDNA$subject_id!="sample8_6" & diversity_long_gDNA$subject_id!="sample8_24"),]

tiff("boxplot_clones_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA_sample8[which(diversity_long_gDNA_sample8$time2==0),], 
            aes(factor(diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==0)]), 
                diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + ylim(0,4000) +
  labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_gDNA_sample8[which(diversity_long_gDNA_sample8$time2==6),], 
            aes(factor(diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==6)]), 
                diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==6)],fill=clin)) + ylim(0,4000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_gDNA_sample8[which(diversity_long_gDNA_sample8$time2==24),], 
            aes(factor(diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==24)]), 
                diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==24)],fill=clin)) + ylim(0,4000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==0)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==0)]))
summary(glm(diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==6)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==6)]))
summary(glm(diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==24)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==24)]))

###If deleted the outlier at time 24 is also signficant?
diversity_long_gDNA_sample6_24<-diversity_long_gDNA[which(diversity_long_gDNA$subject_id!="sample6_24"),]
summary(glm(diversity_long_gDNA_sample6_24$clones_gDNA[which(diversity_long_gDNA_sample6_24$time2==6)] ~ diversity_long_gDNA_sample6_24$clin[which(diversity_long_gDNA_sample6_24$time2==6)]))

tiff("boxplot_entropy_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA_sample8[which(diversity_long_gDNA_sample8$time2==0),], 
            aes(factor(diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==0)]), 
                diversity_long_gDNA_sample8$entropy[which(diversity_long_gDNA_sample8$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + ylim(6,11) +
  labs(title="time 0",x="Clin", y = "entropy")
p2 = ggplot(diversity_long_gDNA_sample8[which(diversity_long_gDNA_sample8$time2==6),], 
            aes(factor(diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==6)]), 
                diversity_long_gDNA_sample8$entropy[which(diversity_long_gDNA_sample8$time2==6)],fill=clin)) + ylim(6,11) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "entropy")
p3 = ggplot(diversity_long_gDNA_sample8[which(diversity_long_gDNA_sample8$time2==24),], 
            aes(factor(diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==24)]), 
                diversity_long_gDNA_sample8$entropy[which(diversity_long_gDNA_sample8$time2==24)],fill=clin)) + ylim(6,11) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "entropy")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA_sample8$entropy[which(diversity_long_gDNA_sample8$time2==0)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==0)]))
summary(glm(diversity_long_gDNA_sample8$entropy[which(diversity_long_gDNA_sample8$time2==6)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==6)]))
summary(glm(diversity_long_gDNA_sample8$entropy[which(diversity_long_gDNA_sample8$time2==24)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==24)]))


###What happened if we adjust for the other covariates?
summary(glm(diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==0)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==0)]+
              diversity_long_gDNA_sample8$Donor.Age[which(diversity_long_gDNA_sample8$time2==0)] +diversity_long_gDNA_sample8$Donor.Source[which(diversity_long_gDNA_sample8$time2==0)] +
              diversity_long_gDNA_sample8$Donor.Gender[which(diversity_long_gDNA_sample8$time2==0)]+diversity_long_gDNA_sample8$recipient.Race[which(diversity_long_gDNA_sample8$time2==0)]  +
              diversity_long_gDNA_sample8$Recipient.Age.when.had.Tx[which(diversity_long_gDNA_sample8$time2==0)]+diversity_long_gDNA_sample8$Recipient.Gender[which(diversity_long_gDNA_sample8$time2==0)]+
              diversity_long_gDNA_sample8$immunosuppression[which(diversity_long_gDNA_sample8$time2==0)]+diversity_long_gDNA_sample8$hla_mismatch[which(diversity_long_gDNA_sample8$time2==0)]))

summary(glm(diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==6)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==6)]+
              diversity_long_gDNA_sample8$Donor.Age[which(diversity_long_gDNA_sample8$time2==6)] +diversity_long_gDNA_sample8$Donor.Source[which(diversity_long_gDNA_sample8$time2==6)] +
              diversity_long_gDNA_sample8$Donor.Gender[which(diversity_long_gDNA_sample8$time2==6)]+diversity_long_gDNA_sample8$recipient.Race[which(diversity_long_gDNA_sample8$time2==6)]  +
              diversity_long_gDNA_sample8$Recipient.Age.when.had.Tx[which(diversity_long_gDNA_sample8$time2==6)]+diversity_long_gDNA_sample8$Recipient.Gender[which(diversity_long_gDNA_sample8$time2==6)]+
              diversity_long_gDNA_sample8$immunosuppression[which(diversity_long_gDNA_sample8$time2==6)]+diversity_long_gDNA_sample8$hla_mismatch[which(diversity_long_gDNA_sample8$time2==6)]))

summary(glm(diversity_long_gDNA_sample8$clones_gDNA[which(diversity_long_gDNA_sample8$time2==24)] ~ diversity_long_gDNA_sample8$clin[which(diversity_long_gDNA_sample8$time2==24)]+
              diversity_long_gDNA_sample8$Donor.Age[which(diversity_long_gDNA_sample8$time2==24)] +diversity_long_gDNA_sample8$Donor.Source[which(diversity_long_gDNA_sample8$time2==24)] +
              diversity_long_gDNA_sample8$Donor.Gender[which(diversity_long_gDNA_sample8$time2==24)]+diversity_long_gDNA_sample8$recipient.Race[which(diversity_long_gDNA_sample8$time2==24)]  +
              diversity_long_gDNA_sample8$Recipient.Age.when.had.Tx[which(diversity_long_gDNA_sample8$time2==24)]+diversity_long_gDNA_sample8$Recipient.Gender[which(diversity_long_gDNA_sample8$time2==24)]+
              diversity_long_gDNA_sample8$immunosuppression[which(diversity_long_gDNA_sample8$time2==24)]+diversity_long_gDNA_sample8$hla_mismatch[which(diversity_long_gDNA_sample8$time2==24)]))


tiff("boxplot_SHM_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + ylim(0.01,0.09) +
  labs(title="time 0",x="Clin", y = "SHM") 
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==6)],fill=clin)) + ylim(0.01,0.09) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "SHM")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==24)],fill=clin)) + ylim(0.01,0.09) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "SHM")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))



########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
fm_null <- lmer(diversity_long_gDNA_sample8$clones_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA_sample8,REML = F)
fm_full <- lmer(diversity_long_gDNA_sample8$clones_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA_sample8,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA_sample8$clones_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()

tiff("Plot_int_clones_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.2) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()


##Diversity
fm_null <- lmer(diversity_long_gDNA_sample8$entropy ~ clin + time + (time | Sample_id),data=diversity_long_gDNA_sample8,REML = F)
fm_full <- lmer(diversity_long_gDNA_sample8$entropy ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA_sample8,REML = F)
anova(fm_full, fm_null) 


tiff("plot_lmer_entropy_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA_sample8$entropy, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5)+
  labs(x = "time (months)",y = "Entropy") 
print(p)
dev.off()

tiff("Plot_int_entropy_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin),size=1.5) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()


##SHM
fm_null <- lmer(diversity_long_gDNA$SHM_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$SHM_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null)

tiff("plot_lmer_SHM_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$SHM_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5)
print(p)
dev.off()
##It is borderline significant

##If we delete the 32 time point
diversity_long_gDNA_32time<-diversity_long_gDNA[which(diversity_long_gDNA$time!=32),]
fm_null <- lmer(diversity_long_gDNA_32time$SHM_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA_32time,REML = F)
fm_full <- lmer(diversity_long_gDNA_32time$SHM_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA_32time,REML = F)
anova(fm_full, fm_null)

tiff("plot_lmer_SHM_gDNA_no32.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA_32time$SHM_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5)
print(p)
dev.off()


################
### AR  gDNA ###
################
diversity_AR<-diversity_reads_clones[which(diversity_reads_clones$clin=="AR" | diversity_reads_clones$clin=="pre-AR"),]
diversity_AR_gDNA<-diversity_AR[which(diversity_AR$reads_gDNA>=100),]
diversity_AR_gDNA$clin<-relevel(diversity_AR_gDNA$clin,ref="pre-AR")

##Number of clones
tiff("plot_summary_clones_AR_gDNA.tiff",h=2500,w=2300,res=300)
ggplot(diversity_AR_gDNA, aes(clin, clones_gDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Clones")
dev.off()

##Entropy
tiff("plot_summary_entropy_AR_gDNA.tiff",h=2500,w=2300,res=300)
ggplot(diversity_AR_gDNA, aes(clin, entropy_gDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Entropy")
dev.off()



##################################
#####Analysis by time and clin ##
##################################
tiff("boxplot_clones_AR_gDNA.tiff",h=2000,w=1800,res=300)
ggplot(diversity_AR_gDNA, aes(clin,clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Clones")
dev.off()
summary(glm(diversity_AR_gDNA$entropy ~ diversity_AR_gDNA$clin))

tiff("boxplot_clones_gDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_gDNA,aes(factor(diversity_AR_gDNA$immunosuppression),
                                  diversity_AR_gDNA$clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Clones")
print(p1)
dev.off()
summary(glm(diversity_AR_gDNA$clones_gDNA ~ diversity_AR_gDNA$clin*diversity_AR_gDNA$immunosuppression))

tiff("boxplot_entropy_gDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_gDNA,aes(factor(diversity_AR_gDNA$immunosuppression),
                                  diversity_AR_gDNA$entropy,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Entropy")
print(p1)
dev.off()
summary(glm(diversity_AR_gDNA$entropy_gDNA ~ diversity_AR_gDNA$clin*diversity_AR_gDNA$immunosuppression))





