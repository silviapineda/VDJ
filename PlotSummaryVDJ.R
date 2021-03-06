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

setwd("/Users/Pinedasans/VDJ/Results/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")


####Plot reads by gDNA 
counts_gDNA <- reads_clones_annot$reads_gDNA[which(reads_clones_annot$reads_gDNA!=0)]
names(counts_gDNA) <-reads_clones_annot$specimen_id[which(reads_clones_annot$reads_gDNA!=0)]
cols = c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod")[reads_clones_annot$clin[which(reads_clones_annot$reads_gDNA!=0)]]

tiff("barplot_reads_by_gDNA.tiff",res=300,w=4000,h=2000)
barplot(counts_gDNA,col=cols,main="Number of Reads gDNA",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend(0, 20000, legend=levels(reads_clones_annot$clin[which(reads_clones_annot$reads_gDNA!=0)]),
       col=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"),pch=15, cex=0.8)
dev.off()

####Plot reads by cDNA 
counts_cDNA <- reads_clones_annot$reads_cDNA[which(reads_clones_annot$reads_cDNA!=0)]
names(counts_cDNA) <-reads_clones_annot$specimen_id[which(reads_clones_annot$reads_cDNA!=0)]
tiff("barplot_reads_by_cDNA.tiff",res=300,w=4000,h=2000)
cols = c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod")[reads_clones_annot$clin[which(reads_clones_annot$reads_cDNA!=0)]]
barplot(counts_cDNA,col=cols,main="Number of Reads cDNA",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend(0, 250000, legend=levels(reads_clones_annot$clin[which(reads_clones_annot$reads_cDNA!=0)]),
       col=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"),pch=15, cex=0.8)
dev.off()

###Plot reads and clones by samples

counts<-table(data_qc$specimen_label,data_qc$amplification_template)
read_count_ighClones<- unique(data_qc[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template")])
counts<-table(read_count_ighClones$specimen_label,read_count_ighClones$amplification_template)

tiff("barplot_reads_clone_by_sample.tiff",res=300,w=4000,h=3000)
par(mfrow = c(2,1))
barplot(t(counts),col=c("darkseagreen","steelblue"),main="Number of Reads",xlab = "Samples", ylab = "Reads",
        legend = colnames(counts),args.legend = list(x = "topleft"),las=2) 

barplot(t(counts),col=c("darkgoldenrod3","brown2"),main="Number of Clones",xlab = "Samples", ylab = "Clones",
        legend = colnames(counts),args.legend = list(x = "topleft"),las=2,ylim=c(0,250000))
dev.off()

###Plot reads and clones by individuals
par(mfrow = c(2,1))
counts<-table(data_qc$sample_id,data_qc$amplification_template)
barplot(t(counts),col=c("darkseagreen","steelblue"),main="Number of Reads",xlab = "Individuals", ylab = "Reads",
        legend = colnames(counts),args.legend = list(x = "topleft"),las=2)

read_count_ighClones<- unique(data_qc[,c("sample_id","V_J_lenghCDR3_Clone_igh","amplification_template")])
counts<-table(read_count_ighClones$sample_id,read_count_ighClones$amplification_template)
barplot(t(counts),col=c("darkgoldenrod3","brown2"),main="Number of Clones",xlab = "Individuals", ylab = "Clones",
        legend = colnames(counts),args.legend = list(x = "topleft"),las=2,ylim=c(0,800000))

## 1. Clones by Cells
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]

##NP
specimen_unique_NP<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$clin=="NP")])
time<-unique(data_qc_gDNA$time[which(data_qc_gDNA$clin=="NP")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_gDNA$specimen_label[which(data_qc_gDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_NP)
  #clones<-table(data_qc_gDNA$V_J_lenghCDR3_Clone_igh,data_qc_gDNA$specimen_label)
  for (i in 1:length(specimen_unique)){
    data_qc_gDNA_specimen<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_gDNA_specimen)[1]>=100){
      clones<-table(data_qc_gDNA_specimen$V_J_lenghCDR3_Clone_igh)
      #write.table(clones,file=paste("NP-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[i],".txt",sep=""))
      tiff(paste("NP-time",time[j],"-ind",data_qc_gDNA_specimen$sample_id[i],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,1400),xlim=c(0,40),col=c("chartreuse4"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
      #ggplot(data.matrix(table(clones)),aes())
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
      plot(data.matrix(table(clones)),ylim=c(0,1400),col=c("firebrick3"),pch=19,xlim=c(0,40),xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
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
      plot(data.matrix(table(clones)),ylim=c(0,1400),col=c("goldenrod"),pch=19,xlim=c(0,40),xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_gDNA_specimen$sample_id[1],",time=",data_qc_gDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}


## 1. Clones by Cells
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

##NP
specimen_unique_NP<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$clin=="NP")])
time<-unique(data_qc_cDNA$time[which(data_qc_cDNA$clin=="NP")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_NP)
  #clones<-table(data_qc_cDNA$V_J_lenghCDR3_Clone_igh,data_qc_cDNA$specimen_label)
  for (i in 1:length(specimen_unique)){
    data_qc_cDNA_specimen<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_cDNA_specimen)[1]>=100){
      clones<-table(data_qc_cDNA_specimen$V_J_lenghCDR3_Clone_igh)
      #write.table(clones,file=paste("NP-time",time[j],"-ind",data_qc_cDNA_specimen$sample_id[i],".txt",sep=""))
      tiff(paste("NP-time",time[j],"-ind",data_qc_cDNA_specimen$sample_id[i],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,22000),xlim=c(0,200),col=c("chartreuse4"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_cDNA_specimen$sample_id[1],",time=",data_qc_cDNA_specimen$time[1],sep=""))
      #ggplot(data.matrix(table(clones)),aes())
      dev.off()
    }
  }
}
##PNR
specimen_unique_PNR<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$clin=="PNR")])
time<-unique(data_qc_cDNA$time[which(data_qc_cDNA$clin=="PNR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_PNR)
  for (i in 1:length(specimen_unique)){
    data_qc_cDNA_specimen<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_cDNA_specimen)[1]>=100){
      clones<-table(data_qc_cDNA_specimen$V_J_lenghCDR3_Clone_igh)
      data.matrix(table(clones))
      tiff(paste("PNR-time",time[j],"-ind",data_qc_cDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,22000),xlim=c(0,200),col=c("dodgerblue3"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_cDNA_specimen$sample_id[1],",time=",data_qc_cDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}
##PR
specimen_unique_PR<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$clin=="PR")])
time<-unique(data_qc_cDNA$time[which(data_qc_cDNA$clin=="PR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_PR)
  for (i in 1:length(specimen_unique)){
    data_qc_cDNA_specimen<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_cDNA_specimen)[1]>=100){
      clones<-table(data_qc_cDNA_specimen$V_J_lenghCDR3_Clone_igh)
      tiff(paste("PR-time",time[j],"-ind",data_qc_cDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,22000),xlim=c(0,200),col=c("darkorange2"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_cDNA_specimen$sample_id[1],",time=",data_qc_cDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}
##AR
specimen_unique_AR<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$clin=="AR")])
time<-unique(data_qc_cDNA$time[which(data_qc_cDNA$clin=="AR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_AR)
  for (i in 1:length(specimen_unique)){
    data_qc_cDNA_specimen<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_cDNA_specimen)[1]>=100){
      clones<-table(data_qc_cDNA_specimen$V_J_lenghCDR3_Clone_igh)
      tiff(paste("AR-time",time[j],"-ind",data_qc_cDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,22000),xlim=c(0,200),col=c("firebrick3"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_cDNA_specimen$sample_id[1],",time=",data_qc_cDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}

##preAR
specimen_unique_preAR<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$clin=="pre-AR")])
time<-unique(data_qc_cDNA$time[which(data_qc_cDNA$clin=="pre-AR")])
time<-time[order(time)]
for (j in 1:length(time)){
  specimen_unique_time<-unique(data_qc_cDNA$specimen_label[which(data_qc_cDNA$time==time[j])])
  specimen_unique<-intersect(specimen_unique_time,specimen_unique_preAR)
  for (i in 1:length(specimen_unique)){
    data_qc_cDNA_specimen<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
    if(dim(data_qc_cDNA_specimen)[1]>=100){
      clones<-table(data_qc_cDNA_specimen$V_J_lenghCDR3_Clone_igh)
      tiff(paste("pre-AR-time",time[j],"-ind",data_qc_cDNA_specimen$sample_id[1],".tiff",sep=""),width = 170, height = 190)
      plot(data.matrix(table(clones)),ylim=c(0,22000),xlim=c(0,200),col=c("goldenrod"),pch=19,xlab = "Number of Cells",ylab = "Number of clones",main = paste("ind=",data_qc_cDNA_specimen$sample_id[1],",time=",data_qc_cDNA_specimen$time[1],sep=""))
      dev.off()
    }
  }
}




###Clones longitudional
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))

##gDNA
reads_clones_annot_Long_gDNA <- reads_clones_annot_Long[which(reads_clones_annot_Long$reads_gDNA>=100),] ## 67 samples (Discard samples with less than 100 reads)
clinLong_gDNA <- clinLong[which(reads_clones_annot_Long$reads_gDNA>=100)] #(Discard samples with less than 100 reads)

p1<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=clones_gDNA, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("ClonesLong gDNA")

##cDNA
reads_clones_annot_Long_cDNA <- reads_clones_annot_Long[which(reads_clones_annot_Long$reads_cDNA>=100),] ## 55 samples (Discard samples with less than 100 reads)
clinLong_cDNA <- clinLong[which(reads_clones_annot_Long$reads_cDNA>=100)] 

p2<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=clones_cDNA, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("ClonesLong cDNA")

##Clones AR
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))

##gDNA
reads_clones_annot_AR_gDNA <- reads_clones_annot_AR[which(reads_clones_annot_AR$reads_gDNA>=100),] ## 8 samples (Discard samples with less than 100 reads)
clinAR_gDNA <- clinAR[which(reads_clones_annot_AR$reads_gDNA>=100)] 

p3<-ggplot(data=reads_clones_annot_AR_gDNA, aes(x=time, y=clones_gDNA, group=subject_id, shape=clinAR_gDNA, color=clinAR_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("ClonesAR gDNA")

##cDNA
reads_clones_annot_AR_cDNA <- reads_clones_annot_AR[which(reads_clones_annot_AR$reads_cDNA>=100),] ## 20 samples (Discard samples with less than 100 reads)
clinAR_cDNA <- clinAR[which(reads_clones_annot_AR$reads_cDNA>=100)] 

p4<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=clones_cDNA, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("ClonesAR cDNA")

ggplot2.multiplot(p1,p2,p3,p4)



################################
##### Demographics for table ###
###############################
table(reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

##Recipient Age
summary(reads_clones_annot$Recipient.Age.when.had.Tx[!duplicated(reads_clones_annot$Sample_id)])
#NP
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="AR")])

fit = lm(reads_clones_annot$Recipient.Age.when.had.Tx[!duplicated(reads_clones_annot$Sample_id)] ~ reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])
anova(fit)



##Donor Age
summary(reads_clones_annot$Donor.Age[!duplicated(reads_clones_annot$Sample_id)])
#NP
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="AR")])

fit = lm(reads_clones_annot$Donor.Age[!duplicated(reads_clones_annot$Sample_id)] ~ reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])
anova(fit)

fit = lm(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)] ~
      reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])
anova(fit)

##immunosuppression
table(reads_clones_annot$immunosuppression[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$immunosuppression[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
      reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])

##Recipient.Gender
table(reads_clones_annot$Recipient.Gender[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$Recipient.Gender[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])


##Donor.Gender
table(reads_clones_annot$Donor.Gender[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$Donor.Gender[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])


##recipient.Race
table(reads_clones_annot$recipient.Race[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$recipient.Race[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])

##Relatedness
table(reads_clones_annot$Donor.Source[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$Donor.Source[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])



###Considering the number of clones gDNA
par(mfrow = c(2, 2))  
boxplot(reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)], 
        main="Time 0", ylab = "Number of clones", col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,5000))
fit = lm((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==0)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)], 
            main="Time 6",ylab = "Number of clones", col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,5000))
fit = lm((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)], ylab = "Number of clones", 
            main="Time 24",col = c("chartreuse4", "dodgerblue3","darkorange2"),ylim=c(0,5000))
fit = lm((reads_clones_annot_Long_gDNA$clones_gDNA[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)])
anova(fit)

#cDNA
boxplot(reads_clones_annot_Long_cDNA$clones_cDNA ~ clinLong_cDNA, main = "Long cDNA" , ylab = "Number of clones",col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_cDNA$clones_cDNA~ clinLong_cDNA))
anova(fit)


####Considering NP as baseline
reads_clones_annot_NP_AR<-reads_clones_annot[which(reads_clones_annot$clin=="NP" | reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
clin_NP_AR<-factor(reads_clones_annot_NP_AR$clin, levels=c("NP", "pre-AR","AR"))

par(mfrow = c(1, 2))  
#gDNA
reads_clones_annot_NP_AR_gDNA<- reads_clones_annot_NP_AR[which(reads_clones_annot_NP_AR$reads_gDNA>=100),] ## 20 samples (Discard samples with less than 100 reads)
clin_NP_AR_gDNA <- clin_NP_AR[which(reads_clones_annot_NP_AR$reads_gDNA>=100)] 
boxplot(reads_clones_annot_NP_AR_gDNA$clones_gDNA ~ clin_NP_AR_gDNA, main = "NP_AR gDNA" , ylab = "Number of clones",col = c("chartreuse4","goldenrod","firebrick3"))
fit = lm((reads_clones_annot_NP_AR_gDNA$clones_gDNA~ clin_NP_AR_gDNA))
anova(fit)


#cDNA
reads_clones_annot_NP_AR_cDNA<- reads_clones_annot_NP_AR[which(reads_clones_annot_NP_AR$reads_cDNA>=100),] ## 20 samples (Discard samples with less than 100 reads)
clin_NP_AR_cDNA <- clin_NP_AR[which(reads_clones_annot_NP_AR$reads_cDNA>=100)] 
boxplot(reads_clones_annot_NP_AR_cDNA$clones_cDNA ~ clin_NP_AR_cDNA, main = "NP_AR cDNA" , ylab = "Number of clones",col = c("chartreuse4","goldenrod","firebrick3"))
fit = lm((reads_clones_annot_NP_AR_cDNA$clones_cDNA~ clin_NP_AR_cDNA))
anova(fit)


