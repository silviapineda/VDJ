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


############################################
####Calculate entropy without downsampling#
###########################################

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


######################################
######## Downsampling the data  ######
######################################
#Downsampling  

##gDNA 1062 reads is the minimum (after QC) for LongData
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
specimen_unique<-unique(data_qc_gDNA$specimen_label)

data_clone_specimen<-NULL
clone_entropy<-matrix(NA,length(specimen_unique),10)
simpson<-matrix(NA,length(specimen_unique),10)
for(j in 1:10){
  print(j)
  for (i in 1:length(specimen_unique)){
    #print(i)
    data_specimen_unique<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(nrow(data_specimen_unique)>=1062){
      data_specimen_unique.down<-data_specimen_unique[sample(nrow(data_specimen_unique),1062),]
      #data_clone_specimen<-rbind(data_clone_specimen,data_specimen_unique.down[,c("specimen_label","igh_clone_id")])
      clone_entropy[i,j]<-entropy(as.numeric(data_specimen_unique.down[,"igh_clone_id"]))
      simpson[i,j]<-simpson(table(as.numeric(data_specimen_unique.down[,"igh_clone_id"])))
      
    }
  }
}

#barplot(table(table(data_clone_specimen$igh_clone_id[which(data_clone_specimen$specimen_label==specimen_unique[6])])))

clone_entropy_mean<-apply(clone_entropy,1,function(x) mean(x,na.rm=T))
clone_entropy_norm<-clone_entropy_mean/max(clone_entropy_mean,na.rm = T)
clonality<-(1-clone_entropy_norm)
names(clonality)<-specimen_unique
simpson_mean<-apply(simpson,1,function(x) mean(x,na.rm=T))
names(simpson_mean)<-specimen_unique


id.sample<-match(names(clonality),reads_clones_annot$specimen_id)
reads_clones_annot[id.sample,"clonality"]<-clonality
reads_clones_annot[id.sample,"simpson"]<-simpson_mean


##Separate Long and AR
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]

clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))

fm1 <- lmer(clonality~1 + (1 | Sample_id),data=reads_clones_annot_Long)
fm2 <- lmer(clonality~ time*clin + (time | Sample_id),data=reads_clones_annot_Long)
# also works with lme objects
anova(fm1, fm2) 


####Boxplot and ANOVA analysis for Longitudinal Data
###gDNA

p1<-ggplot(data=reads_clones_annot_Long, aes(x=time, y=clonality, group=subject_id, shape=clinLong, color=clinLong)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality Long gDNA")

# 
# reads_clones_annot_Long.noNA<-reads_clones_annot_Long[which(is.na(reads_clones_annot_Long$clonality)==F),]
# data<-data.frame(reads_clones_annot_Long.noNA[which(reads_clones_annot_Long.noNA$time==0 |
#                                                       reads_clones_annot_Long.noNA$time==6 |
#                                                       reads_clones_annot_Long.noNA$time==24),c("Sample_id","clin","time","clonality")])
# ggplot(data, aes(x=time, y=clonality, fill=Sample_id)) + 
#    geom_area()

par(mfrow = c(2, 2))  
boxplot(reads_clones_annot_Long$clonality[which(reads_clones_annot_Long$time==0)] ~ clinLong[which(reads_clones_annot_Long$time==0)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main="Time 0")
boxplot(reads_clones_annot_Long$clonality[which(reads_clones_annot_Long$time==6)] ~ clinLong[which(reads_clones_annot_Long$time==6)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main = "Time 6")
boxplot(reads_clones_annot_Long$clonality[which(reads_clones_annot_Long$time==24)] ~ clinLong[which(reads_clones_annot_Long$time==24)] , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"),main="Time 24")


# boxplot(reads_clones_annot_Long$clonality ~ clinLong , ylab = "Clonality", 
#         col = c("chartreuse4", "dodgerblue3","darkorange2"))
# fit = lm(reads_clones_annot_Long$clonality ~ clinLong)
# anova(fit)

p2<-ggplot(data=reads_clones_annot_AR, aes(x=time, y=clonality, group=subject_id, shape=clinAR, color=clinAR)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Clonality AR gDNA")

boxplot(reads_clones_annot_AR$clonality ~ clinAR , ylab = "Clonality", 
        col = c("goldenrod","firebrick3"))
fit = lm(reads_clones_annot_AR$clonality ~ clinAR)
anova(fit)



##cDNA 5124 reads is the minimum (after QC) for ARdata
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
specimen_unique<-unique(data_qc_cDNA$specimen_label)

data_specimen_unique.down<-NULL
clone_entropy_cDNA<-matrix(NA,length(specimen_unique),10)
simpson_cDNA<-matrix(NA,length(specimen_unique),10)
for(j in 1:10){
  print(j)
  for (i in 1:length(specimen_unique)){
    data_specimen_unique<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
    if(nrow(data_specimen_unique)>=5124){
      data_specimen_unique.down<-data_specimen_unique[sample(nrow(data_specimen_unique),5124),]
      clone_entropy_cDNA[i,j]<-entropy(as.numeric(data_specimen_unique.down[,"igh_clone_id"]))
      simpson_cDNA[i,j]<-simpson(table(as.numeric(data_specimen_unique.down[,"igh_clone_id"])))
      
    }
  }
}

clone_entropy_mean_cDNA<-apply(clone_entropy_cDNA,1,function(x) mean(x,na.rm=T))
clone_entropy_norm_cDNA<-clone_entropy_mean_cDNA/max(clone_entropy_mean_cDNA,na.rm = T)
clonality_cDNA<-(1-clone_entropy_norm_cDNA)
names(clonality_cDNA)<-specimen_unique

id.sample<-match(names(clonality_cDNA),reads_clones_annot$specimen_id)
reads_clones_annot[id.sample,"clonality_cDNA"]<-clonality_cDNA


##Separate Long and AR
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]

clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))

reads_clones_annot_AR_baseline<-reads_clones_annot[which(reads_clones_annot$clin=="AR"| reads_clones_annot$clin=="pre-AR" | reads_clones_annot$clin=="NP"),]
clinAR_baseline<-factor(reads_clones_annot_AR_baseline$clin, levels=c("pre-AR","AR","NP"))
####Boxplot and ANOVA analysis for Longitudinal Data
###gDNA

###Fitthing a longitudinal model
fm1 <- lmer(clonality_cDNA~time + (time | Sample_id),data=reads_clones_annot_Long)
fm2 <- lmer(clonality_cDNA~ time*clin + (time | Sample_id),data=reads_clones_annot_Long)
# also works with lme objects
anova(fm1, fm2)

###Fitthing a longitudinal model
fm1 <- lmer(clonality_cDNA~time + (time | Sample_id),data=reads_clones_annot_AR)
fm2 <- lmer(clonality_cDNA~ clin + (1 | Sample_id),data=reads_clones_annot_AR)
# also works with lme objects
anova(fm1, fm2)

p1.cDNA<-ggplot(data=reads_clones_annot_Long, aes(x=time, y=clonality_cDNA, group=subject_id, shape=clinLong, color=clinLong)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Clonality Long cDNA")


boxplot(reads_clones_annot_Long$clonality_cDNA ~ clinLong , ylab = "Clonality", 
        col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm(reads_clones_annot_Long$clonality_cDNA ~ clinLong)
anova(fit)

p2.cDNA<-ggplot(data=reads_clones_annot_AR, aes(x=time, y=clonality_cDNA, group=subject_id, shape=clinAR, color=clinAR)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("Clonality AR cDNA")

boxplot(reads_clones_annot_AR$clonality_cDNA ~ clinAR , ylab = "Clonality", 
        col = c("goldenrod","firebrick3"))
fit = lm(reads_clones_annot_AR$clonality_cDNA ~ clinAR)
anova(fit)

p3<-ggplot(data=reads_clones_annot_AR_baseline, aes(x=time, y=clonality_cDNA, group=subject_id, shape=clinAR_baseline, color=clinAR_baseline)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3","chartreuse4")) +
  ggtitle("Clonality AR baseline NP cDNA")

