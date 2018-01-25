
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


#########
## Calculate demographics
########

reads_clones_annot_reads100<-reads_clones_annot[which(reads_clones_annot$total_reads>100),]
reads_clones_annot_reads100_long<-reads_clones_annot_reads100[which(reads_clones_annot_reads100$clin!="pre-AR" & reads_clones_annot_reads100$clin!="AR"),]
reads_clones_annot_reads100_AR<-reads_clones_annot_reads100[which(reads_clones_annot_reads100$clin=="pre-AR" | reads_clones_annot_reads100$clin=="AR"),]

##reads_clones_annot_subjects<-rbind(reads_clones_annot_reads100_long[!duplicated(reads_clones_annot_reads100_long$Sample_id),],
#                                   reads_clones_annot_reads100_AR[!duplicated(reads_clones_annot_reads100_AR$Sample_id),])

#Only longitudinal data
reads_clones_annot_subjects<-rbind(reads_clones_annot_reads100_long[!duplicated(reads_clones_annot_reads100_long$Sample_id),])

dim(reads_clones_annot_subjects)
table(reads_clones_annot_subjects$clin)

##CADI score
summary(reads_clones_annot_subjects$cadi)
summary(reads_clones_annot_subjects$cadi[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$cadi[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$cadi[which(reads_clones_annot_subjects$clin=="PR")])
summary(reads_clones_annot_subjects$cadi[which(reads_clones_annot_subjects$clin=="pre-AR" | reads_clones_annot_subjects$clin=="AR")])

fit = lm(reads_clones_annot_subjects$cadi ~ reads_clones_annot_subjects$clin)
anova(fit)


##Donor Age
summary(reads_clones_annot_subjects$Donor.Age)
summary(reads_clones_annot_subjects$Donor.Age[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$Donor.Age[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$Donor.Age[which(reads_clones_annot_subjects$clin=="PR")])
summary(reads_clones_annot_subjects$Donor.Age[which(reads_clones_annot_subjects$clin=="pre-AR" | reads_clones_annot_subjects$clin=="AR")])

fit = lm(reads_clones_annot_subjects$Donor.Age ~ reads_clones_annot_subjects$clin)
anova(fit)

##Recipient Age
summary(reads_clones_annot_subjects$Recipient.Age.when.had.Tx)
summary(reads_clones_annot_subjects$Recipient.Age.when.had.Tx[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$Recipient.Age.when.had.Tx[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$Recipient.Age.when.had.Tx[which(reads_clones_annot_subjects$clin=="PR")])
summary(reads_clones_annot_subjects$Recipient.Age.when.had.Tx[which(reads_clones_annot_subjects$clin=="pre-AR" | reads_clones_annot_subjects$clin=="AR")])

fit = lm(reads_clones_annot_subjects$Recipient.Age.when.had.Tx ~ reads_clones_annot_subjects$clin)
anova(fit)

##Donor Gender
summary(reads_clones_annot_subjects$Donor.Gender)
summary(reads_clones_annot_subjects$Donor.Gender[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$Donor.Gender[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$Donor.Gender[which(reads_clones_annot_subjects$clin=="PR")])

chisq.test(rbind(c(7,5,2),c(3,5,5)))

##Recipient Gender
summary(reads_clones_annot_subjects$Recipient.Gender)
summary(reads_clones_annot_subjects$Recipient.Gender[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$Recipient.Gender[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$Recipient.Gender[which(reads_clones_annot_subjects$clin=="PR")])

chisq.test(rbind(c(5,6,5),c(5,4,2)))

##Donor Source
summary(reads_clones_annot_subjects$Donor.Source)
summary(reads_clones_annot_subjects$Donor.Source[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$Donor.Source[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$Donor.Source[which(reads_clones_annot_subjects$clin=="PR")])

chisq.test(rbind(c(6,4,3),c(4,6,4)))

##Reciepient race
summary(reads_clones_annot_subjects$recipient.Race)
summary(reads_clones_annot_subjects$recipient.Race[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$recipient.Race[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$recipient.Race[which(reads_clones_annot_subjects$clin=="PR")])

chisq.test(rbind(c(2,6,4),c(4,1,1),c(0,0,1),c(4,3,1)))

##Immunosuprrsion
summary(reads_clones_annot_subjects$immunosuppression)
summary(reads_clones_annot_subjects$immunosuppression[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$immunosuppression[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$immunosuppression[which(reads_clones_annot_subjects$clin=="PR")])

chisq.test(rbind(c(6,4,3),c(4,6,4)))

##HLA
summary(reads_clones_annot_subjects$hla_mismatch)
summary(reads_clones_annot_subjects$hla_mismatch[which(reads_clones_annot_subjects$clin=="NP")])
summary(reads_clones_annot_subjects$hla_mismatch[which(reads_clones_annot_subjects$clin=="PNR")])
summary(reads_clones_annot_subjects$hla_mismatch[which(reads_clones_annot_subjects$clin=="PR")])

fit = lm(reads_clones_annot_subjects$hla_mismatch ~ reads_clones_annot_subjects$clin)
anova(fit)


##table time by outcome
table(reads_clones_annot_reads100$clin[which(reads_clones_annot_reads100$reads_gDNA!=0)],reads_clones_annot_reads100$time[which(reads_clones_annot_reads100$reads_gDNA!=0)])
table(reads_clones_annot_reads100$clin[which(reads_clones_annot_reads100$reads_cDNA!=0)],reads_clones_annot_reads100$time[which(reads_clones_annot_reads100$reads_cDNA!=0)])


########################################
####Calculate repertoire diversity ####
#######################################

##########
## gDNA ##
##########
data_gDNA<-data_merge[which(data_merge$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin=="NP" | data_gDNA$clin=="PNR" | data_gDNA$clin=="PR"),]

###########################
## 2. Diversity measures###
###########################
############
## gDNA ### 
###########
specimen_unique<-unique(data_gDNA$specimen_label)
entropy<-rep(NA,length(specimen_unique))
simpson<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_gDNA[which(data_gDNA$specimen_label==specimen_unique[i]),]
  clones_specimen<-data_specimen_unique[,"V_J_lenghCDR3_CloneId"]
  #To write file to run with Recon
  #write.delim(data.frame(table(table(clones_specimen))),file=paste("clones_",specimen_unique[i],".txt",sep=""),sep="\t",col.names=F)
  fi<-as.numeric(table(clones_specimen))/length(clones_specimen)
  hi<-fi*log2(fi)
  entropy[i]=-sum(hi) #entropy(table(clones_specimen)) returns the same result but by default is the natural logarithm (log)
  simpson[i]=sum(fi*fi) 
}
entropy_norm<-entropy/max(entropy,na.rm = T)
clonality<-(1-entropy_norm)
names(clonality)<-specimen_unique
diversity<-cbind(clonality,entropy,simpson)

write.csv(diversity,"/Users/Pinedasans/VDJ/Data/diversity_AllClones_gDNA.csv")


########################################
####  Statistical Analysis  ############
#######################################

diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity_AllClones_gDNA.csv",header=T)

reads_clones_annot_gDNA<-reads_clones_annot[which(reads_clones_annot$gDNA=="gDNA"),]
id<-match(reads_clones_annot_gDNA$specimen_id,diversity[,1])
diversity_reads_clones<-cbind(reads_clones_annot_gDNA,diversity[id,-1])

#############
### gDNA ####
#############

diversity_gDNA<-diversity_reads_clones[which(diversity_reads_clones$reads_gDNA>=100),]

###Longitudinal Data
diversity_long_gDNA<-diversity_gDNA[which(diversity_gDNA$clin=="NP" | diversity_gDNA$clin=="PNR" | diversity_gDNA$clin=="PR"),]
diversity_long_gDNA$clin<-factor(diversity_long_gDNA$clin, levels=c("NP", "PNR", "PR"))
table(diversity_long_gDNA$clin[which(is.na(diversity_long_gDNA$reads_gDNA)==F)],diversity_long_gDNA$time[which(is.na(diversity_long_gDNA$reads_gDNA)==F)])

ggplot(data=diversity_long_gDNA, aes(x=time, y=clones_PCR_A, group=Sample_id, shape=clin, color=clin)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal gDNA")


##################################
#####Analysis by time and clin ##
##################################

diversity_long_gDNA$time2<-replace(diversity_long_gDNA$time,diversity_long_gDNA$time==12,6)

###Barplots
#clones
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
tiff("barplot_clones_time_gDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==0)],
        cex.names=0.8,las=2,ylim = c(0,4000),ylab = c("Clones"))
legend(0, 4000, legend=levels(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==6)],
        cex.names=0.8,las=2,ylim = c(0,4000),ylab = c("Clones"))
barplot(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==24)],
        cex.names=0.8,las=2,ylim = c(0,4000),ylab = c("Clones"))
dev.off()

#reads
tiff("barplot_reads_time_gDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_gDNA$reads_gDNA[which(diversity_long_gDNA$time2==0)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==0)],
        cex.names=0.8,las=2,ylim = c(0,20000),ylab = c("Reads"))
legend(0, 20000, legend=levels(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_gDNA$reads_gDNA[which(diversity_long_gDNA$time2==6)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==6)],
        cex.names=0.8,las=2,ylim = c(0,20000),ylab = c("Reads"))
barplot(diversity_long_gDNA$reads_gDNA[which(diversity_long_gDNA$time2==24)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==24)],
        cex.names=0.8,las=2,ylim = c(0,20000),ylab = c("Reads"))
dev.off()


#SHM
diversity_long_gDNA$SHM_gDNA_byClones<-diversity_long_gDNA$SHM_gDNA/diversity_long_gDNA$clones_gDNA
tiff("barplot_SHM_byClones_gDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==0)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==0)],
        cex.names=0.8,las=2,ylim = c(0,0.2),ylab = c("SHM"))
legend(0, .2, legend=levels(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==6)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==6)],
        cex.names=0.8,las=2,ylim = c(0,0.2),ylab = c("SHM"))
barplot(diversity_long_gDNA$SHM_gDNA_byClones[which(diversity_long_gDNA$time2==24)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==24)],
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
  labs(x = "time (months)",y = "Entropy (gDNA)") 
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





