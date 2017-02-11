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
library(gplots)
require(graphics)
library(RColorBrewer)
library("bcRep")
library(ggplot2)
library(devtools)
library(easyGgplot2)

setwd("/Users/Pinedasans/Documents/VDJ/")
load("/Users/Pinedasans/Data/VDJ/VDJ.Rdata")

###################################
### Study with No Downsampling ####
##################################

##### 1.Longitudinal data
setwd("/Users/Pinedasans/Documents/VDJ/NoDownSamplingStudy/")
## 1.1plot reads
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))

##gDNA
reads_clones_annot_Long_gDNA <- reads_clones_annot_Long[which(reads_clones_annot_Long$gDNA_reads>=100),] ## 67 samples (Discard samples with less than 100 reads)
clinLong_gDNA <- clinLong[which(reads_clones_annot_Long$gDNA_reads>=100)] #(Discard samples with less than 100 reads)

p1.gDNA<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=gDNA_reads, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Total-Reads Long gDNA")

##cDNA
reads_clones_annot_Long_cDNA <- reads_clones_annot_Long[which(reads_clones_annot_Long$cDNA_reads>=100),] ## 67 samples (Discard samples with less than 100 reads)
clinLong_cDNA <- clinLong[which(reads_clones_annot_Long$cDNA_reads>=100)] #(Discard samples with less than 100 reads)

p1.cDNA<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=cDNA_reads, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Total-Reads Long cDNA")

## 1.2 plot clones

#gDNA
p2.gDNA<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=clones_gDNA, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Total-Clones Long gDNA")

#cDNA
p2.cDNA<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=clones_cDNA, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Total-Clones Long cDNA")

ggplot2.multiplot(p1.gDNA,p2.gDNA,p1.cDNA,p2.cDNA)

## 1.3 boxplot - ANOVA
#gDNA
par(mfrow = c(2, 2))  
boxplot(reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)], 
        main="Time 0 gDNA", ylab = "Number of clones", col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,5000))
fit = lm((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==0)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)], 
        main="Time 6 gDNA",ylab = "Number of clones", col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,5000))
fit = lm((reads_clones_annot_Long_gDNA$clones_cDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)], ylab = "Number of clones", 
        main="Time 24 gDNA",col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,5000))
fit = lm((reads_clones_annot_Long_gDNA$clones_cDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)])
anova(fit)

#cDNA
par(mfrow = c(1, 2))  
boxplot(reads_clones_annot_Long_cDNA$clones_cDNA[which(reads_clones_annot_Long_cDNA$time==6)] ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==6)], 
        main="Time 6 cDNA", ylab = "Number of clones", col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,70000))
fit = lm((reads_clones_annot_Long_cDNA$clones_cDNA[which(reads_clones_annot_Long_cDNA$time==6)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_cDNA$clones_cDNA[which(reads_clones_annot_Long_cDNA$time==24)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==24)], ylab = "Number of clones", 
        main="Time 24 cDNA",col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,70000))
fit = lm((reads_clones_annot_Long_cDNA$clones_cDNA[which(reads_clones_annot_Long_cDNA$time==24)]) ~ clinLong_cDNA[which(reads_clones_annot_Long_cDNA$time==24)])
anova(fit)

##### 2.AR

## 2.1 plot reads
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))

##gDNA
reads_clones_annot_AR_gDNA <- reads_clones_annot_AR[which(reads_clones_annot_AR$gDNA_reads>=100),] ## 10 samples (Discard samples with less than 100 reads)
clinAR_gDNA <- clinAR[which(reads_clones_annot_AR$gDNA_reads>=100)] 

p1.gDNA<-ggplot(data=reads_clones_annot_AR_gDNA, aes(x=time, y=gDNA_reads, group=subject_id, shape=clinAR_gDNA, color=clinAR_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Total-Reads AR gDNA")

##cDNA
reads_clones_annot_AR_cDNA <- reads_clones_annot_AR[which(reads_clones_annot_AR$cDNA_reads>=100),] ## 10 samples (Discard samples with less than 100 reads)
clinAR_cDNA <- clinAR[which(reads_clones_annot_AR$cDNA_reads>=100)] 

p1.cDNA<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=cDNA_reads, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Total-Reads AR cDNA")

## 2.2 plot clones

#gDNA
p2.gDNA<-ggplot(data=reads_clones_annot_AR_gDNA, aes(x=time, y=clones_gDNA, group=subject_id, shape=clinAR_gDNA, color=clinAR_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Total-Clones AR gDNA")

#cDNA
p2.cDNA<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=clones_cDNA, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Total-Clones AR cDNA")

ggplot2.multiplot(p1.gDNA,p2.gDNA,p1.cDNA,p2.cDNA)

## 2.3 boxplot with NP as baseline- ANOVA
reads_clones_annot_NP_AR<-reads_clones_annot[which(reads_clones_annot$clin=="NP" | reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
clin_NP_AR<-factor(reads_clones_annot_NP_AR$clin, levels=c("NP", "pre-AR","AR"))

par(mfrow = c(1, 2))  
#gDNA
reads_clones_annot_NP_AR_gDNA<- reads_clones_annot_NP_AR[which(reads_clones_annot_NP_AR$gDNA_reads>=100),] ## 20 samples (Discard samples with less than 100 reads)
clin_NP_AR_gDNA <- clin_NP_AR[which(reads_clones_annot_NP_AR$gDNA_reads>=100)] 
boxplot(reads_clones_annot_NP_AR_gDNA$clones_gDNA ~ clin_NP_AR_gDNA, main = "NP_AR gDNA" , ylab = "Number of clones",col = c("chartreuse4","goldenrod","firebrick3"))
fit = lm((reads_clones_annot_NP_AR_gDNA$clones_gDNA~ clin_NP_AR_gDNA))
anova(fit)


#cDNA
reads_clones_annot_NP_AR_cDNA<- reads_clones_annot_NP_AR[which(reads_clones_annot_NP_AR$cDNA_reads>=100),] ## 20 samples (Discard samples with less than 100 reads)
clin_NP_AR_cDNA <- clin_NP_AR[which(reads_clones_annot_NP_AR$cDNA_reads>=100)] 
boxplot(reads_clones_annot_NP_AR_cDNA$clones_cDNA ~ clin_NP_AR_cDNA, main = "NP_AR cDNA" , ylab = "Number of clones",col = c("chartreuse4","goldenrod","firebrick3"))
fit = lm((reads_clones_annot_NP_AR_cDNA$clones_cDNA~ clin_NP_AR_cDNA))
anova(fit)
