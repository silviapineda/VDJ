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

setwd("/Users/Pinedasans/VDJ/Results/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")

#########
## Calculate demographics
########

reads_clones_annot_reads100<-reads_clones_annot[which(reads_clones_annot$total_reads>100),]
reads_clones_annot_reads100_long<-reads_clones_annot_reads100[which(reads_clones_annot_reads100$clin!="pre-AR" & reads_clones_annot_reads100$clin!="AR"),]
reads_clones_annot_reads100_AR<-reads_clones_annot_reads100[which(reads_clones_annot_reads100$clin=="pre-AR" | reads_clones_annot_reads100$clin=="AR"),]

reads_clones_annot_subjects<-rbind(reads_clones_annot_reads100_long[!duplicated(reads_clones_annot_reads100_long$Sample_id),],
                                   reads_clones_annot_reads100_AR[!duplicated(reads_clones_annot_reads100_AR$Sample_id),])

dim(reads_clones_annot_subjects)
table(reads_clones_annot_subjects$clin)

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


##table time by outcome
table(reads_clones_annot_reads100$clin[which(reads_clones_annot_reads100$reads_gDNA!=0)],reads_clones_annot_reads100$time[which(reads_clones_annot_reads100$reads_gDNA!=0)])
table(reads_clones_annot_reads100$clin[which(reads_clones_annot_reads100$reads_cDNA!=0)],reads_clones_annot_reads100$time[which(reads_clones_annot_reads100$reads_cDNA!=0)])



####################################################
### Compare same individuals from gDNA and cDNA ###
###################################################


########################################
####Calculate repertoire diversity ####
#######################################

##########
## gDNA ##
##########
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]

##########
## cDNA ##
##########
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]


###########################
## 2. Diversity measures###
###########################
############
## gDNA ### 
###########
specimen_unique<-unique(data_qc_gDNA$specimen_label)
entropy<-rep(NA,length(specimen_unique))
simpson<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
  clones_specimen<-data_specimen_unique[,"V_J_lenghCDR3_Clone_igh"]
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

write.csv(diversity,"diversity_gDNA.csv")

######Read the output from RECON with the estimates Hill numbers
results<-read.table("/Users/Pinedasans/VDJ/Results/clones_D_number_table.txt", header = T,sep = "\t")
entropy<-log(results$obs_1.0D)
entropy_estimates<-log(results$est_1.0D)
simpson<-1/results$obs_2.0D
simpson_estimates<-1/results$est_2.0D

entropy_norm<-entropy_estimates/max(entropy_estimates,na.rm = T)
clonality<-(1-entropy_norm)

clones<-results$est_0.0D
diversity_estimates<-cbind(clonality,entropy_estimates,simpson_estimates,clones)
rownames(diversity_estimates)<-specimen_unique[order(specimen_unique)]
write.csv(diversity_estimates,"/Users/Pinedasans/VDJ/Results/diversity_estimates_gDNA.csv")


##########
## cDNA ##
##########

## 1. Clones by Cells
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

specimen_unique<-unique(data_qc_cDNA$specimen_label)
entropy<-rep(NA,length(specimen_unique))
simpson<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  #print(i)
  data_specimen_unique<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
  clones_specimen<-data_specimen_unique[,"V_J_lenghCDR3_Clone_igh"]
  ##To apply Recon
  #write.delim(data.frame(table(table(clones_specimen))),file=paste("clones_cDNA",specimen_unique[i],".txt",sep=""),sep="\t",col.names=F)
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

######Read the output from RECON with the estimates Hill numbers
results<-read.table("/Users/Pinedasans/VDJ/Results/test_D_number_table_cDNA.txt", header = T,sep = "\t")
entropy<-log(results$obs_1.0D)
entropy_estimates<-log(results$est_1.0D)
simpson<-1/results$obs_2.0D
simpson_estimates<-1/results$est_2.0D

entropy_norm<-entropy_estimates/max(entropy_estimates,na.rm = T)
clonality<-(1-entropy_norm)

clones<-results$est_0.0D
diversity_estimates<-cbind(clonality,entropy_estimates,simpson_estimates,clones)
rownames(diversity_estimates)<-specimen_unique[order(specimen_unique)]
write.csv(diversity_estimates,"/Users/Pinedasans/VDJ/Results/diversity_estimates_cDNA.csv")


########################################
####  Statistical Analysis  ############
#######################################
diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity.csv",header=T)


#############
### gDNA ####
#############
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
diversity_gDNA_qc<-diversity[which(diversity$reads_gDNA>=100),]

cor(diversity_gDNA_qc$clones_gDNA,diversity_gDNA_qc$clones_est_gDNA,use = "complete.obs")
tiff("plot_cor_clones_gDNA.tiff",res=300,w=1500,h=1500)
qplot(clones_gDNA, clones_est_gDNA, data = diversity_gDNA_qc, colour = clin) +
      scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity_gDNA_qc$entropy_gDNA,diversity_gDNA_qc$entropy_est_gDNA,use = "complete.obs")
tiff("plot_cor_entropy_gDNA.tiff",res=300,w=1500,h=1500)
qplot(entropy_gDNA, entropy_est_gDNA, data = diversity_gDNA_qc, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()


cor(diversity_gDNA_qc$simpson_gDNA,diversity_gDNA_qc$simpson_est_gDNA,use = "complete.obs")
tiff("plot_cor_simpson_gDNA.tiff",res=300,w=1500,h=1500)
qplot(simpson_gDNA, simpson_est_gDNA, data = diversity_gDNA_qc, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity_gDNA_qc$clonality_gDNA,diversity_gDNA_qc$clonality_est_gDNA,use = "complete.obs")
tiff("plot_cor_clonality_gDNA.tiff",res=300,w=1500,h=1500)
qplot(clonality_gDNA, clonality_est_gDNA, data = diversity_gDNA_qc, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

###Longitudinal Data
diversity_long_gDNA<-diversity_gDNA_qc[which(diversity_gDNA_qc$clin=="NP" | diversity_gDNA_qc$clin=="PNR" | diversity_gDNA_qc$clin=="PR"),]
diversity_long_gDNA$clin<-factor(diversity_long_gDNA$clin, levels=c("NP", "PNR", "PR"))
table(diversity_long_gDNA$clin[which(is.na(diversity_long_gDNA$reads_gDNA)==F)],diversity_long_gDNA$time[which(is.na(diversity_long_gDNA$reads_gDNA)==F)])

p1<-ggplot(data=diversity_long_gDNA, aes(x=time, y=clones_gDNA, group=Sample_id, shape=clin, color=clin)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal gDNA")

##Number of clones
g1<-ggplot(diversity_long_gDNA[1:24,], aes(time, clones_gDNA)) + scale_y_continuous(limit = c(0,5000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_gDNA[25:48,], aes(time, clones_gDNA)) + geom_point() + scale_y_continuous(limit = c(0,5000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_gDNA[49:69,], aes(time, clones_gDNA)) + geom_point() + scale_y_continuous(limit = c(0,5000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")

tiff("plot_summary_clones_gDNA.tiff",h=2500,w=2300,res=300)
multiplot(g1, g2, g3, rows=3)
dev.off()

##Entropy
g1<-ggplot(diversity_long_gDNA[1:24,], aes(time, entropy_gDNA)) + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Entropy")
g2<-ggplot(diversity_long_gDNA[25:48,], aes(time, entropy_gDNA)) + geom_point() + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Entropy")
g3<-ggplot(diversity_long_gDNA[49:69,], aes(time, entropy_gDNA)) + geom_point() + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Entropy")
tiff("plot_summary_entropy_gDNA.tiff",h=2500,w=2300,res=300)
multiplot(g1, g2, g3, rows=3)
dev.off()

##Simpson
g1<-ggplot(diversity_long_gDNA[1:24,], aes(time, simpson_gDNA)) + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_gDNA[25:48,], aes(time, simpson_gDNA)) + geom_point() + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_gDNA[49:69,], aes(time, simpson_gDNA)) + geom_point() + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
multiplot(g1, g2, g3, rows=3)


##correlation 
cor(diversity_long_gDNA$reads_gDNA,diversity_long_gDNA$clonality_gDNA)
qplot(clones_gDNA, clonality_gDNA, data = diversity_long_gDNA, colour = clin) +
  scale_colour_manual(values=c("chartreuse4","dodgerblue3","darkorange2"))



##################################
#####Analysis by time and clin ##
##################################
tiff("boxplot_clones_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))

summary(glm(diversity_long_gDNA$clones_gDNA ~ diversity_long_gDNA$clin))


tiff("boxplot_entropy_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "entropy")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "entropy")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "entropy")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))

tiff("boxplot_clonality_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$clonality_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "clonality")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$clonality_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "clonality")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$clonality_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "clonality")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$clonality_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))

tiff("boxplot_simpson_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "simpson")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "simpson")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "simpson")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))


tiff("boxplot_SHM_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$SHM_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "SHM")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$SHM_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "SHM")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$SHM_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "SHM")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$SHM_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$SHM_gDNA[which(diversity_long_gDNA$time2==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$SHM_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))

tiff("boxplot_CDR3length_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "cdr3_length")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "cdr3_length")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "cdr3_length")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))



########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
fm_null <- lmer(diversity_long_gDNA$clones_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$clones_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

#ggplot(fm_full, aes(time, diversity_long_gDNA$clones_gDNA, color=clin)) +
#  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
#  stat_summary(fun.data=mean_se, geom="pointrange") +
# stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
#  geom_line(aes(y = predict(fm_full))) 

tiff("plot_lmer_clones_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5)
  #geom_line(aes(eff_df$time,eff_df$fit,group=clin,colour=clin),size=1.5) 
print(p)
dev.off()

tiff("Plot_int1_clones_gDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_clones_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.2) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Diversity
fm_null <- lmer(diversity_long_gDNA$entropy_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$entropy_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_entropy_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$entropy_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5)
print(p)
dev.off()

tiff("Plot_int1_entropy_gDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_entropy_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin),size=1.5) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Diversity
fm_null <- lmer(diversity_long_gDNA$ ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$entropy_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

##SHM
fm_null <- lmer(diversity_long_gDNA$SHM_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$SHM_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 
tiff("Plot_int2_SHM_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin),size=1.5) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

#################################
###Considering other variables##
################################

###All together 
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]+
          diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==0)]+diversity_long_gDNA$Recipient.Gender[which(diversity_long_gDNA$time2==0)]+
  diversity_long_gDNA$Donor.Gender[which(diversity_long_gDNA$time2==0)]+diversity_long_gDNA$GenderMismatch[which(diversity_long_gDNA$time2==0)]+
  diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==0)]+diversity_long_gDNA$Donor.Age[which(diversity_long_gDNA$time2==0)]+
  diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==0)]+diversity_long_gDNA$recipient.Race[which(diversity_long_gDNA$time2==0)]+
    diversity_long_gDNA$hla_mismatch[which(diversity_long_gDNA$time2==0)]+diversity_long_gDNA$cadi[which(diversity_long_gDNA$time2==0)]))
##PR at time 0 is still significant after adjustment of all the variables

summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]+
              diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==6)]+diversity_long_gDNA$Recipient.Gender[which(diversity_long_gDNA$time2==6)]+
              diversity_long_gDNA$Donor.Gender[which(diversity_long_gDNA$time2==6)]+diversity_long_gDNA$GenderMismatch[which(diversity_long_gDNA$time2==6)]+
              diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==6)]+diversity_long_gDNA$Donor.Age[which(diversity_long_gDNA$time2==6)]+
              diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==6)]+diversity_long_gDNA$recipient.Race[which(diversity_long_gDNA$time2==6)]+
              diversity_long_gDNA$hla_mismatch[which(diversity_long_gDNA$time2==6)]+diversity_long_gDNA$cadi[which(diversity_long_gDNA$time2==6)]))

summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]+
              diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==24)]+diversity_long_gDNA$Recipient.Gender[which(diversity_long_gDNA$time2==24)]+
              diversity_long_gDNA$Donor.Gender[which(diversity_long_gDNA$time2==24)]+diversity_long_gDNA$GenderMismatch[which(diversity_long_gDNA$time2==24)]+
              diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==24)]+diversity_long_gDNA$Donor.Age[which(diversity_long_gDNA$time2==24)]+
              diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==24)]+diversity_long_gDNA$recipient.Race[which(diversity_long_gDNA$time2==24)]+
              diversity_long_gDNA$hla_mismatch[which(diversity_long_gDNA$time2==24)]+diversity_long_gDNA$cadi[which(diversity_long_gDNA$time2==24)]))


##Considering all variables as confounding factors in the clin interaction
fm_null <- lmer(clones_gDNA ~ clin + time + immunosuppression + Recipient.Gender + Donor.Gender + GenderMismatch + Recipient.Age.when.had.Tx + Donor.Age
                + Donor.Source + recipient.Race + hla_mismatch + cadi + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_gDNA ~  clin*time  + immunosuppression + Recipient.Gender + Donor.Gender + GenderMismatch + Recipient.Age.when.had.Tx + Donor.Age
                + Donor.Source + recipient.Race + hla_mismatch + cadi + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_adj_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##############
#### 1. Immuno
#############
COLOR=brewer.pal(3,"Set1")

##adjust only for immuno
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]+
              diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==0)]))
###Not a confounding factor

###Asscoiation with immuno
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)] ~ +
              diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)] ~ +
              diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)] ~ +
              diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==24)]))

summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)] ~ +
              diversity_long_gDNA$hla_mismatch[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)] ~ +
              diversity_long_gDNA$hla_mismatch[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)] ~ +
              diversity_long_gDNA$hla_mismatch[which(diversity_long_gDNA$time2==24)]))
###Plot 
tiff("boxplot_entropy_gDNA_Immuno.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)],fill=immunosuppression)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 0",x="Immunosuppression", y = "Entropy")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)],fill=immunosuppression)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 6",x="Immunosuppression", y = "Entropy")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)],fill=immunosuppression)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 24",x="Immunosuppression", y = "Entropy")
multiplot(p1,p2,p3)
dev.off()


###Association with clin stratified by immunosupression
tiff("boxplot_CDR3_gDNA_Clin_by_immuno.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Immunosuppression", y = "mean_CDR3_length")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Immunosuppression", y = "mean_CDR3_length")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==24)]), 
              diversity_long_gDNA$mean_CDR3_length_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Immunosuppression", y = "mean_CDR3_length")
multiplot(p1,p2,p3)
dev.off()

tiff("boxplot_entropy_gDNA_Clin_by_immuno.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Immunosuppression", y = "Entropy")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Immunosuppression", y = "Entropy")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Immunosuppression", y = "Entropy")
multiplot(p1,p2,p3)
dev.off()


####Longitudinal model for immunosuppression by time
fm_null <- lmer(clones_gDNA ~  time +immunosuppression + (time | Sample_id),
                data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_gDNA ~  time*immunosuppression + (time | Sample_id) ,
                data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_immuno_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("immunosuppression", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=immunosuppression)) +
  geom_line(aes(time,fit,group=immunosuppression,colour=immunosuppression),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=COLOR)
dev.off()

fm_null <- lmer(entropy_gDNA ~  time +immunosuppression + (time | Sample_id),
                data=diversity_long_gDNA,REML = F)
fm_full <- lmer(entropy_gDNA ~  time*immunosuppression + (time | Sample_id) ,
                data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_entropy_immuno_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("immunosuppression", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=immunosuppression)) +
  geom_line(aes(time,fit,group=immunosuppression,colour=immunosuppression),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=COLOR)
dev.off()

###Adjust the association with clin by immunosuppression
fm_null <- lmer(entropy_gDNA ~ clin + time + immunosuppression + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(entropy_gDNA ~  clin*time + immunosuppression + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_entropy_clin_by_immuno_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Interaction cllin*immuno*time
fm_null <- lmer(clones_gDNA ~ clin + time + immunosuppression +time*clin + clin*immunosuppression +time*immunosuppression + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_gDNA ~  clin*time*immunosuppression + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int_clones_clin_immuno_time_gDNA.tiff",h=1600,w=2700,res=300)
clones.effects <- allEffects(fm_full)
plot(clones.effects, ylab="Clones", rug=FALSE)
dev.off()

#CDR3-length
fm_null <- lmer(mean_CDR3_length_gDNA ~ clin + time + immunosuppression +time*clin + clin*immunosuppression +time*immunosuppression + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(mean_CDR3_length_gDNA ~  clin*time*immunosuppression + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int_CDR3_clin_immuno_time_gDNA.tiff",h=1600,w=2700,res=300)
clones.effects <- allEffects(fm_full)
plot(clones.effects, ylab="Clones", rug=FALSE)
dev.off()


########
## 2. Donor.Source 
########
COLOR=brewer.pal(3,"Set1")

##adjust only for donor source
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]+
              diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==0)]))
###Not a confounding factor

###Asscoiation with donor source
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)] ~ +
              diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)] ~ +
              diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)] ~ +
              diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==24)]))


###Plot 
tiff("boxplot_clones_gDNA_Donor.Source.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 0",x="Donor.Source", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 6",x="Donor.Source", y = "Clones")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 24",x="Donor.Source", y = "Clones")
multiplot(p1,p2,p3)
dev.off()

tiff("boxplot_entropy_gDNA_Donor.Source.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 0",x="Donor.Source", y = "entropy")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 6",x="Donor.Source", y = "entropy")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 24",x="Donor.Source", y = "entropy")
multiplot(p1,p2,p3)
dev.off()

###Association with clin stratified by immunosupression
tiff("boxplot_clones_gDNA_Clin_by_DonorSource.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Donor.Source", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Donor.Source", y = "Clones")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Donor.Source", y = "Clones")
multiplot(p1,p2,p3)
dev.off()

tiff("boxplot_entropy_gDNA_Clin_by_Donor.Source.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==0)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Donor.Source", y = "Entropy")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Donor.Source", y = "Entropy")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Donor.Source", y = "Entropy")
multiplot(p1,p2,p3)
dev.off()


####Longitudinal model for immunosuppression by time
fm_null <- lmer(entropy_gDNA ~  time +Donor.Source + (time | Sample_id),
                data=diversity_long_gDNA,REML = F)
fm_full <- lmer(entropy_gDNA ~  time*Donor.Source + (time | Sample_id) ,
                data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_entropy_Donor.Source_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("Donor.Source", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=Donor.Source)) +
  geom_line(aes(time,fit,group=Donor.Source,colour=Donor.Source),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=COLOR)
dev.off()


###Adjust the association with clin by immunosuppression
fm_null <- lmer(entropy_gDNA ~ clin + time + Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(entropy_gDNA ~  clin*time + Donor.Source+ (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_entropy_clin_by_donorSource_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Interaction cllin*immuno*time
fm_null <- lmer(clones_gDNA ~ clin*time + clin*Donor.Source+time*Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_gDNA ~  clin*time*Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int_entropy_clin_donorSource_time_gDNA.tiff",h=1600,w=2700,res=300)
clones.effects <- allEffects(fm_full)
plot(clones.effects, ylab="Clones", rug=FALSE)
dev.off()



########
## 3. Recipient Age 
########
##adjust only for donor source
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]+
              diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==0)]))
###Not a confounding factor

###Asscoiation with donor source
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)] ~ +
              diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)] ~ +
              diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)] ~ +
              diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==24)]))

###Association with RecAge
tiff("plot_clones_gDNA_Rec_Age.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),],aes(diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==0)], 
            diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)]))+geom_point() + geom_smooth(method='lm') + labs(title="time 0",x="Rec.Age", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),],aes(diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==6)], 
           diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)]))+geom_point() + geom_smooth(method='lm') + labs(title="time 6",x="Rec.Age", y = "Clones")
p3 =ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),],aes(diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time2==24)], 
           diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)]))+geom_point() + geom_smooth(method='lm') + labs(title="time 24",x="Rec.Age", y = "Clones")

multiplot(p1,p2,p3)
dev.off()



####Longitudinal model for RecAge by time
fm_null <- lmer(clones_gDNA ~  time +Recipient.Age.when.had.Tx + (time | Sample_id),
                data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_gDNA ~  time*Recipient.Age.when.had.Tx + (time | Sample_id) ,
                data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 



###Adjust the association with clin by Donor.Source
fm_null <- lmer(clones_gDNA ~ clin + time + Recipient.Age.when.had.Tx + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_gDNA ~  clin*time + Recipient.Age.when.had.Tx + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 


tiff("Plot_int2_clones_clin_by_recAge_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()


##Interaction cllin*immuno*time
fm_null <- lmer(entropy_gDNA ~ clin*time + clin*Recipient.Age.when.had.Tx+time*Recipient.Age.when.had.Tx + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(entropy_gDNA ~  clin*time*Recipient.Age.when.had.Tx + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 


###Adjust the association with clin by Donor.Source and immuno
fm_null <- lmer(clones_gDNA ~ clin + time + immunosuppression + Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_gDNA ~  clin*time + immunosuppression + Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_clin_by_immuno_DonorSource_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()




################
### AR  gDNA ###
################
diversity_AR<-diversity[which(diversity$clin=="AR" | diversity$clin=="pre-AR"),]
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
summary(glm(diversity_AR_gDNA$entropy_gDNA ~ diversity_AR_gDNA$clin))

tiff("boxplot_clones_gDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_gDNA,aes(factor(diversity_AR_gDNA$immunosuppression),
                                  diversity_AR_gDNA$clones_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Clones")
print(p1)
dev.off()
summary(glm(diversity_AR_gDNA$clones_gDNA ~ diversity_AR_gDNA$clin*diversity_AR_gDNA$immunosuppression))

tiff("boxplot_entropy_gDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_gDNA,aes(factor(diversity_AR_gDNA$immunosuppression),
                                  diversity_AR_gDNA$entropy_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Entropy")
print(p1)
dev.off()
summary(glm(diversity_AR_gDNA$entropy_gDNA ~ diversity_AR_gDNA$clin*diversity_AR_gDNA$immunosuppression))


tiff("boxplot_SHM_gDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_gDNA,aes(factor(diversity_AR_gDNA$immunosuppression),
                                  diversity_AR_gDNA$SHM_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "SHM")
print(p1)
dev.off()


tiff("boxplot_CDR3_gDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_gDNA,aes(factor(diversity_AR_gDNA$immunosuppression),
                                  diversity_AR_gDNA$mean_CDR3_length_gDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "mean_CDR3_length")
print(p1)
dev.off()





