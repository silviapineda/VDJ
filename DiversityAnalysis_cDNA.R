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

####################################################
### Compare same individuals from gDNA and cDNA ###
###################################################


########################################
####Calculate repertoire diversity ####
#######################################

##########
## cDNA ##
##########
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]



########################################
####  Statistical Analysis  ############
#######################################
diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity.csv",header=T)


#############
### cDNA ####
#############
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
diversity_cDNA_qc<-diversity[which(diversity$reads_cDNA>=100),]


cor(diversity_cDNA_qc$clones_cDNA,diversity_cDNA_qc$clones_est_cDNA,use = "complete.obs")
tiff("plot_cor_clones_cDNA.tiff",res=300,w=1500,h=1500)
qplot(clones_cDNA, clones_est_cDNA, data = diversity_cDNA_qc, colour = clin) +
      scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity_cDNA_qc$entropy_cDNA,diversity_cDNA_qc$entropy_est_cDNA,use = "complete.obs")
tiff("plot_cor_entropy_cDNA.tiff",res=300,w=1500,h=1500)
qplot(entropy_cDNA, entropy_est_cDNA, data = diversity_cDNA_qc, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()


cor(diversity_cDNA_qc$simpson_cDNA,diversity_cDNA_qc$simpson_est_cDNA,use = "complete.obs")
tiff("plot_cor_simpson_cDNA.tiff",res=300,w=1500,h=1500)
qplot(simpson_cDNA, simpson_est_cDNA, data = diversity_cDNA_qc, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity_cDNA_qc$clonality_cDNA,diversity_cDNA_qc$clonality_est_cDNA,use = "complete.obs")
tiff("plot_cor_clonality_cDNA.tiff",res=300,w=1500,h=1500)
qplot(clonality_cDNA, clonality_est_cDNA, data = diversity_cDNA_qc, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

###Longitudinal Data
diversity_long_cDNA<-diversity_cDNA_qc[which(diversity_cDNA_qc$clin=="NP" | diversity_cDNA_qc$clin=="PNR" | diversity_cDNA_qc$clin=="PR"),]
diversity_long_cDNA$clin<-factor(diversity_long_cDNA$clin, levels=c("NP", "PNR", "PR"))
table(diversity_long_cDNA$clin[which(is.na(diversity_long_cDNA$reads_cDNA)==F)],diversity_long_cDNA$time[which(is.na(diversity_long_cDNA$reads_cDNA)==F)])

diversity_long_cDNA_noOutlier<-diversity_long_cDNA[which(diversity_long_cDNA$Sample_id!="300003"),]

p1<-ggplot(data=diversity_long_cDNA, aes(x=time, y=clones_cDNA, group=Sample_id, shape=clin, color=clin)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal cDNA")

##Number of clones
g1<-ggplot(diversity_long_cDNA[1:18,], aes(time, clones_cDNA)) + scale_y_continuous(limit = c(0,80000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_cDNA[19:38,], aes(time, clones_cDNA)) + geom_point() + scale_y_continuous(limit = c(0,80000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_cDNA[39:55,], aes(time, clones_cDNA)) + geom_point() +  scale_y_continuous(limit = c(0,80000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")

tiff("plot_summary_clones_cDNA.tiff",h=2500,w=2300,res=300)
multiplot(g1, g2, g3, rows=3)
dev.off()

##Entropy
g1<-ggplot(diversity_long_cDNA[1:18,], aes(time, entropy_cDNA)) + scale_y_continuous(limit = c(7,15)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_cDNA[19:38,], aes(time, entropy_cDNA)) + geom_point() + scale_y_continuous(limit = c(7,15)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_cDNA[39:55,], aes(time, entropy_cDNA)) + geom_point() + scale_y_continuous(limit = c(7,15)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
tiff("plot_summary_entropy_cDNA.tiff",h=2500,w=2300,res=300)
multiplot(g1, g2, g3, rows=3)
dev.off()

##Simpson
g1<-ggplot(diversity_long_cDNA[1:18,], aes(time, simpson_cDNA)) + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_cDNA[19:38,], aes(time, simpson_cDNA)) + geom_point() + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_cDNA[39:55,], aes(time, simpson_cDNA)) + geom_point() + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
tiff("plot_summary_simpson_cDNA.tiff",h=2500,w=2300,res=300)
multiplot(g1, g2, g3, rows=3)
dev.off()

##################################
#####Analysis by time and clin ##
##################################
tiff("boxplot_clones_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==6)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]))
summary(glm(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==24)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]))

tiff("boxplot_entropy_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "entropy")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "entropy")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "entropy")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]))
summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]))


tiff("boxplot_SHM_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$SHM_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "SHM")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$SHM_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "SHM")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$SHM_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "SHM")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_cDNA$SHM[which(diversity_long_cDNA$time2==24)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]))

tiff("boxplot_CDR3length_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$mean_CDR3_length_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "cdr3_length")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$mean_CDR3_length_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "cdr3_length")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$mean_CDR3_length_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "cdr3_length")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_cDNA$mean_CDR3_length[which(diversity_long_cDNA$time2==24)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]))



########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
fm_null <- lmer(diversity_long_cDNA_noOutlier$clones_cDNA ~ clin + time + (time | Sample_id),data=diversity_long_cDNA_noOutlier,REML = F)
fm_full <- lmer(diversity_long_cDNA_noOutlier$clones_cDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_cDNA_noOutlier,REML = F)
anova(fm_full, fm_null) 

#ggplot(fm_full, aes(time, diversity_long_cDNA$clones_cDNA, color=clin)) +
#  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
#  stat_summary(fun.data=mean_se, geom="pointrange") +
# stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
#  geom_line(aes(y = predict(fm_full))) 

tiff("plot_lmer_clones_cDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_cDNA_noOutlier$clones_cDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5)
  #geom_line(aes(eff_df$time,eff_df$fit,group=clin,colour=clin),size=1.5) 
print(p)
dev.off()

tiff("Plot_int1_clones_cDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_SHM_cDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.2) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Diversity
fm_null <- lmer(diversity_long_cDNA_noOutlier$entropy_cDNA ~ clin + time + (time | Sample_id),data=diversity_long_cDNA_noOutlier,REML = F)
fm_full <- lmer(diversity_long_cDNA_noOutlier$entropy_cDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_cDNA_noOutlier,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_entropy_cDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_cDNA$entropy_cDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5)
print(p)
dev.off()

tiff("Plot_int1_entropy_cDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_entropy_cDNA.tiff",h=1000,w=1500,res=300)
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
summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]+
              diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==6)]+diversity_long_cDNA$Recipient.Gender[which(diversity_long_cDNA$time2==6)]+
              diversity_long_cDNA$Donor.Gender[which(diversity_long_cDNA$time2==6)]+diversity_long_cDNA$GenderMismatch[which(diversity_long_cDNA$time2==6)]+
              diversity_long_cDNA$Recipient.Age.when.had.Tx[which(diversity_long_cDNA$time2==6)]+diversity_long_cDNA$Donor.Age[which(diversity_long_cDNA$time2==6)]+
              diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==6)]+diversity_long_cDNA$recipient.Race[which(diversity_long_cDNA$time2==6)]+
              diversity_long_cDNA$hla_mismatch[which(diversity_long_cDNA$time2==6)]+diversity_long_cDNA$cadi[which(diversity_long_cDNA$time2==6)]))

summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]+
              diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==24)]+diversity_long_cDNA$Recipient.Gender[which(diversity_long_cDNA$time2==24)]+
              diversity_long_cDNA$Donor.Gender[which(diversity_long_cDNA$time2==24)]+diversity_long_cDNA$GenderMismatch[which(diversity_long_cDNA$time2==24)]+
              diversity_long_cDNA$Recipient.Age.when.had.Tx[which(diversity_long_cDNA$time2==24)]+diversity_long_cDNA$Donor.Age[which(diversity_long_cDNA$time2==24)]+
              diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==24)]+diversity_long_cDNA$recipient.Race[which(diversity_long_cDNA$time2==24)]+
              diversity_long_cDNA$hla_mismatch[which(diversity_long_cDNA$time2==24)]+diversity_long_cDNA$cadi[which(diversity_long_cDNA$time2==24)]))


##Considering all variables as confounding factors in the clin interaction
fm_null <- lmer(clones_cDNA ~ clin + time + immunosuppression + Recipient.Gender + Donor.Gender + GenderMismatch + Recipient.Age.when.had.Tx + Donor.Age
                + Donor.Source + recipient.Race + hla_mismatch + cadi + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(clones_cDNA ~  clin*time  + immunosuppression + Recipient.Gender + Donor.Gender + GenderMismatch + Recipient.Age.when.had.Tx + Donor.Age
                + Donor.Source + recipient.Race + hla_mismatch + cadi + (time | Sample_id),data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_adj_cDNA.tiff",h=1000,w=1500,res=300)
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

###Asscoiation with immuno
summary(glm(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==6)] ~ +
              diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==6)]))
summary(glm(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==24)] ~ +
              diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==24)]))

summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)] ~ +
              diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==6)]))
summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)] ~ +
              diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==24)]))

###Plot 
tiff("boxplot_entropy_cDNA_Immuno.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==0)],fill=immunosuppression)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 0",x="Immunosuppression", y = "Entropy")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)],fill=immunosuppression)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 6",x="Immunosuppression", y = "Entropy")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)],fill=immunosuppression)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 24",x="Immunosuppression", y = "Entropy")
multiplot(p1,p2,p3)
dev.off()


###Association with clin stratified by immunosupression
tiff("boxplot_clones_cDNA_Clin_by_immuno.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Immunosuppression", y = "Clones")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Immunosuppression", y = "Clones")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Immunosuppression", y = "Clones")
multiplot(p1,p2,p3)
dev.off()

tiff("boxplot_entropy_cDNA_Clin_by_immuno.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Immunosuppression", y = "Entropy")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Immunosuppression", y = "Entropy")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$immunosuppression[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Immunosuppression", y = "Entropy")
multiplot(p1,p2,p3)
dev.off()


####Longitudinal model for immunosuppression by time
fm_null <- lmer(clones_cDNA ~  time +immunosuppression + (time | Sample_id),
                data=diversity_long_cDNA,REML = F)
fm_full <- lmer(clones_cDNA ~  time*immunosuppression + (time | Sample_id) ,
                data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_immuno_cDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("immunosuppression", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=immunosuppression)) +
  geom_line(aes(time,fit,group=immunosuppression,colour=immunosuppression),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=COLOR)
dev.off()


###Adjust the association with clin by immunosuppression
fm_null <- lmer(clones_cDNA ~ clin + time + immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(clones_cDNA ~  clin*time + immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_clin_by_immuno_cDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Interaction cllin*immuno*time
fm_null <- lmer(SHM_cDNA ~ clin + time + immunosuppression +time*clin + clin*immunosuppression +time*immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(SHM_cDNA ~  clin*time*immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int_SHM_clin_immuno_time_cDNA.tiff",h=1600,w=2700,res=300)
clones.effects <- allEffects(fm_full)
plot(clones.effects, ylab="SHM", rug=FALSE)
dev.off()

##entropy
fm_null <- lmer(entropy_cDNA ~ clin + time + immunosuppression +time*clin + clin*immunosuppression +time*immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(entropy_cDNA ~  clin*time*immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int_entropy_clin_immuno_time_cDNA.tiff",h=1600,w=2700,res=300)
clones.effects <- allEffects(fm_full)
plot(clones.effects, ylab="Clones", rug=FALSE)
dev.off()

##entropy
fm_null <- lmer(mean_CDR3_length_cDNA ~ clin + time + immunosuppression +time*clin + clin*immunosuppression +time*immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(mean_CDR3_length_cDNA ~  clin*time*immunosuppression + (time | Sample_id),data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int_CDR3_clin_immuno_time_cDNA.tiff",h=1600,w=2700,res=300)
clones.effects <- allEffects(fm_full)
plot(clones.effects, ylab="mean_CDR3_length", rug=FALSE)
dev.off()

########
## 2. Donor.Source 
########
COLOR=brewer.pal(3,"Set1")


###Asscoiation with donor source
summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)] ~ +
              diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==6)]))
summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)] ~ +
              diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==24)]))


###Plot 
tiff("boxplot_clones_cDNA_Donor.Source.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==0)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 0",x="Donor.Source", y = "Clones")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==6)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 6",x="Donor.Source", y = "Clones")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==24)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 24",x="Donor.Source", y = "Clones")
multiplot(p1,p2,p3)
dev.off()

tiff("boxplot_entropy_cDNA_Donor.Source.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==0)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 0",x="Donor.Source", y = "entropy")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 6",x="Donor.Source", y = "entropy")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)],fill=Donor.Source)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 24",x="Donor.Source", y = "entropy")
multiplot(p1,p2,p3)
dev.off()

###Association with clin stratified by immunosupression
tiff("boxplot_clones_cDNA_Clin_by_DonorSource.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Donor.Source", y = "Clones")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Donor.Source", y = "Clones")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Donor.Source", y = "Clones")
multiplot(p1,p2,p3)
dev.off()

tiff("boxplot_entropy_cDNA_Clin_by_Donor.Source.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 0",x="Donor.Source", y = "Entropy")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 6",x="Donor.Source", y = "Entropy")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$Donor.Source[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))  + labs(title="time 24",x="Donor.Source", y = "Entropy")
multiplot(p1,p2,p3)
dev.off()


####Longitudinal model for immunosuppression by time
fm_null <- lmer(clones_cDNA ~  time +Donor.Source + (time | Sample_id),
                data=diversity_long_cDNA,REML = F)
fm_full <- lmer(clones_cDNA ~  time*Donor.Source + (time | Sample_id) ,
                data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_Donor.Source_cDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("Donor.Source", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=Donor.Source)) +
  geom_line(aes(time,fit,group=Donor.Source,colour=Donor.Source),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=COLOR)
dev.off()


###Adjust the association with clin by immunosuppression
fm_null <- lmer(clones_cDNA ~ clin + time + Donor.Source + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(clones_cDNA ~  clin*time + Donor.Source+ (time | Sample_id),data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int2_clones_clin_by_donorSource_cDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin),size=1.5)+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Interaction cllin*immuno*time
fm_null <- lmer(entropy_cDNA ~ clin*time + clin*Donor.Source+time*Donor.Source + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(entropy_cDNA ~  clin*time*Donor.Source + (time | Sample_id),data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("Plot_int_entropy_clin_donorSource_time_cDNA.tiff",h=1600,w=2700,res=300)
clones.effects <- allEffects(fm_full)
plot(clones.effects, ylab="Clones", rug=FALSE)
dev.off()



################
### AR  cDNA ###
################
diversity_AR<-diversity[which(diversity$clin=="AR" | diversity$clin=="pre-AR"),]
diversity_AR_cDNA<-diversity_AR[which(diversity_AR$reads_cDNA>=100),]
diversity_AR_cDNA$clin<-relevel(diversity_AR_cDNA$clin,ref="pre-AR")

##Number of clones
tiff("plot_summary_clones_AR_cDNA.tiff",h=2500,w=2300,res=300)
ggplot(diversity_AR_cDNA, aes(clin, clones_cDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Clones")
dev.off()

##Entropy
tiff("plot_summary_entropy_AR_cDNA.tiff",h=2500,w=2300,res=300)
ggplot(diversity_AR_cDNA, aes(clin, entropy_cDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Entropy")
dev.off()



##################################
#####Analysis by time and clin ##
##################################
tiff("boxplot_clones_AR_cDNA.tiff",h=2000,w=1800,res=300)
ggplot(diversity_AR_cDNA, aes(clin,clones_cDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Clones")
dev.off()
summary(glm(diversity_AR_cDNA$clones_cDNA ~ diversity_AR_cDNA$clin))

tiff("boxplot_entropy_AR_cDNA.tiff",h=2000,w=1800,res=300)
ggplot(diversity_AR_cDNA, aes(clin,entropy_cDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Clones")
dev.off()
summary(glm(diversity_AR_cDNA$entropy_cDNA ~ diversity_AR_cDNA$clin))

tiff("boxplot_clones_cDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_cDNA,aes(factor(diversity_AR_cDNA$immunosuppression),
                                  diversity_AR_cDNA$SHM_cDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Clones")
print(p1)
dev.off()

summary(glm(diversity_AR_cDNA$SHM_cDNA ~ diversity_AR_cDNA$clin*diversity_AR_cDNA$immunosuppression))


tiff("boxplot_entropy_cDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_cDNA,aes(factor(diversity_AR_cDNA$immunosuppression),
                                  diversity_AR_cDNA$entropy_cDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Entropy")
print(p1)
dev.off()

summary(glm(diversity_AR_cDNA$entropy_cDNA ~ diversity_AR_cDNA$clin*diversity_AR_cDNA$immunosuppression))

tiff("boxplot_SHM_cDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_cDNA,aes(factor(diversity_AR_cDNA$immunosuppression),
                                  diversity_AR_cDNA$SHM_cDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "SHM")
print(p1)
dev.off()


tiff("boxplot_CDR3_cDNA_Clin_by_immuno_AR.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_AR_cDNA,aes(factor(diversity_AR_cDNA$immunosuppression),
                                  diversity_AR_cDNA$mean_CDR3_length_cDNA,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "mean_CDR3_length")
print(p1)
dev.off()
summary(glm(diversity_AR_cDNA$mean_CDR3_length_cDNA ~ diversity_AR_cDNA$clin*diversity_AR_cDNA$immunosuppression))


