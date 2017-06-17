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
diversity<-read.table("/Users/Pinedasans/VDJ/Data/diversity.txt",sep="\t",header=T)


#############
### gDNA ####
#############
cor(diversity$clones_gDNA,diversity$clones_est_gDNA,use = "complete.obs")
tiff("plot_cor_clones_gDNA.tiff",res=300,w=1500,h=1500)
qplot(clones_gDNA, clones_est_gDNA, data = diversity, colour = clin) +
      scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity$entropy_gDNA,diversity$entropy_est_gDNA,use = "complete.obs")
tiff("plot_cor_entropy_gDNA.tiff",res=300,w=1500,h=1500)
qplot(entropy_gDNA, entropy_est_gDNA, data = diversity, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity$simpson_gDNA,diversity$simpson_est_gDNA,use = "complete.obs")
tiff("plot_cor_simpson_gDNA.tiff",res=300,w=1500,h=1500)
qplot(simpson_gDNA, simpson_est_gDNA, data = diversity, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()


###Longitudinal Data
diversity_long<-diversity[which(diversity$clin=="NP" | diversity$clin=="PNR" | diversity$clin=="PR"),]
diversity_long_gDNA<-diversity_long[which(diversity_long$reads_gDNA>=100),]
clinLong<-factor(diversity_long_gDNA$clin, levels=c("NP", "PNR", "PR"))
table(diversity_long_qc$clin[which(is.na(diversity_long_qc$reads_gDNA)==F)],diversity_long_qc$time[which(is.na(diversity_long_qc$reads_gDNA)==F)])

p1<-ggplot(data=diversity_long_gDNA, aes(x=time, y=clones_est_gDNA, group=Sample_id, shape=clinLong, color=clinLong)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal gDNA")

##Number of clones
g1<-ggplot(diversity_long_gDNA[1:21,], aes(time, clones_gDNA)) + scale_y_continuous(limit = c(0,5000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_gDNA[22:48,], aes(time, clones_gDNA)) + geom_point() + scale_y_continuous(limit = c(0,5000)) + 
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
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_gDNA[25:48,], aes(time, entropy_gDNA)) + geom_point() + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_gDNA[49:69,], aes(time, entropy_gDNA)) + geom_point() + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
multiplot(g1, g2, g3, rows=3)

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

##################################
#####Analysis by time and clin ##
##################################
tiff("plot_violin_clones_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==0)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==0)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==6)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==6)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==24)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==24)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 24",x="Clin", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time==0)]))

tiff("plot_violin_entropy_est_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==0)]), 
                diversity_long_gDNA$entropy_est_gDNA[which(diversity_long_gDNA$time==0)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 0",x="Clin", y = "Entropy")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==6)]), 
                diversity_long_gDNA$entropy_est_gDNA[which(diversity_long_gDNA$time==6)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 6",x="Clin", y = "Entropy")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==24)]), 
                diversity_long_gDNA$entropy_est_gDNA[which(diversity_long_gDNA$time==24)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 24",x="Clin", y = "Entropy")

multiplot(p1,p2,p3)
dev.off()

summary(glm(diversity_long_gDNA$entropy_gDNA[which(diversity_long_gDNA$time==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time==6)]))

tiff("plot_violin_simpson_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==0)]), 
                diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time==0)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 0",x="Clin", y = "Entropy")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==6)]), 
                diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time==6)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 6",x="Clin", y = "Entropy")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time==24)]), 
                diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time==24)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 24",x="Clin", y = "Entropy")

multiplot(p1,p2,p3)
dev.off()

summary(glm(diversity_long_gDNA$simpson_gDNA[which(diversity_long_gDNA$time==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time==0)]))

########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
fm_null <- lmer(diversity_long_gDNA$clones_est_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$clones_est_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

#ggplot(fm_full, aes(time, diversity_long_gDNA$clones_gDNA, color=clin)) +
#  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
#  stat_summary(fun.data=mean_se, geom="pointrange") +
# stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
#  geom_line(aes(y = predict(fm_full))) 

tiff("plot_lmer_clones_est_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_est_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int1_clones_est_gDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_clones_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Diversity
fm_null <- lmer(diversity_long_gDNA$entropy_est_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$entropy_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_entropy_est_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int1_entropy_est_gDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_entropy_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

#################################
###Considering other variables##
################################
summary(glm(clones_gDNA~Donor.Source,data=diversity_long_gDNA))
COLOR=brewer.pal(3,"Set1")


########
## 1. Immunosuppression 
########
tiff("boxplot_immuno_clones_gDNA.tif",h=2000,w=2000,res=300)
boxplot(diversity_long_gDNA$clones_gDNA~diversity_long_gDNA$immunosuppression,col=COLOR,main="clones by immunosupression")
dev.off()
summary(glm(clones_gDNA~immunosuppression,data=diversity_long_gDNA))
summary(glm(clones_gDNA~clin,data=diversity_long_gDNA))

###Association with immunosupression
tiff("plot_violin_clones_gDNA_immunosupression.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==0),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time==0)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==0)],fill=immunosuppression)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=COLOR) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 0",x="immunosuppression", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==6),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time==6)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==6)],fill=immunosuppression)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=COLOR) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 6",x="immunosuppression", y = "Clones")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==24),], 
            aes(factor(diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time==24)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==24)],fill=immunosuppression)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=COLOR) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 24",x="immunosuppression", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==24)] ~ diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time==24)]))

##Association with clin adjusted by immuno
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time==0)]
            + diversity_long_gDNA$immunosuppression[which(diversity_long_gDNA$time==0)]))

###Association with clin stratified by immunosupression
diversity_long_gDNA_steFree<-diversity_long_gDNA[which(diversity_long_gDNA$immunosuppression=="Steroid-free"),]
tiff("plot_violin_clones_gDNA_clin_steFree.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA_steFree[which(diversity_long_gDNA_steFree$time==0),], 
            aes(factor(diversity_long_gDNA_steFree$clin[which(diversity_long_gDNA_steFree$time==0)]), 
                diversity_long_gDNA_steFree$clones_gDNA[which(diversity_long_gDNA_steFree$time==0)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_gDNA_steFree[which(diversity_long_gDNA_steFree$time==6),], 
            aes(factor(diversity_long_gDNA_steFree$clin[which(diversity_long_gDNA_steFree$time==6)]), 
                diversity_long_gDNA_steFree$clones_gDNA[which(diversity_long_gDNA_steFree$time==6)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_gDNA_steFree[which(diversity_long_gDNA_steFree$time==24),], 
            aes(factor(diversity_long_gDNA_steFree$clin[which(diversity_long_gDNA_steFree$time==24)]), 
                diversity_long_gDNA_steFree$clones_gDNA[which(diversity_long_gDNA_steFree$time==24)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 24",x="Clin", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA_steFree$clones_gDNA[which(diversity_long_gDNA_steFree$time==24)] ~ diversity_long_gDNA_steFree$clin[which(diversity_long_gDNA_steFree$time==24)]))

diversity_long_gDNA_steBased<-diversity_long_gDNA[which(diversity_long_gDNA$immunosuppression=="Steroid-based"),]
tiff("plot_violin_clones_gDNA_clin_steBase.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA_steBased[which(diversity_long_gDNA_steBased$time==0),], 
            aes(factor(diversity_long_gDNA_steBased$clin[which(diversity_long_gDNA_steBased$time==0)]), 
                diversity_long_gDNA_steBased$clones_gDNA[which(diversity_long_gDNA_steBased$time==0)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_gDNA_steBased[which(diversity_long_gDNA_steBased$time==6),], 
            aes(factor(diversity_long_gDNA_steBased$clin[which(diversity_long_gDNA_steBased$time==6)]), 
                diversity_long_gDNA_steBased$clones_gDNA[which(diversity_long_gDNA_steBased$time==6)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_gDNA_steBased[which(diversity_long_gDNA_steBased$time==24),], 
            aes(factor(diversity_long_gDNA_steBased$clin[which(diversity_long_gDNA_steBased$time==24)]), 
                diversity_long_gDNA_steBased$clones_gDNA[which(diversity_long_gDNA_steBased$time==24)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 24",x="Clin", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA_steBased$clones_gDNA[which(diversity_long_gDNA_steBased$time==24)] ~ diversity_long_gDNA_steBased$clin[which(diversity_long_gDNA_steBased$time==24)]))

####Longitudinal model for immunosuppression by time
fm_null <- lmer(clones_est_gDNA ~  time +immunosuppression + (time | Sample_id),
                data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_est_gDNA ~  time*immunosuppression + (time | Sample_id) ,
                data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_immuno_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_gDNA, colour = immunosuppression)) +
  scale_colour_manual(values=COLOR) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int2_clones_immuno_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("immunosuppression", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=immunosuppression)) +
  geom_line(aes(time,fit,group=immunosuppression,colour=immunosuppression))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=COLOR)
dev.off()


###Adjust the association with clin by immunosuppression
fm_null <- lmer(clones_est_gDNA ~ clin + time + Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_est_gDNA ~  clin*time + Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_clin_by_immuno_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int2_clones_clin_by_immuno_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

########
## 2. Donor.Source 
########
tiff("boxplot_Donor.Source_clones_gDNA.tif",h=2000,w=2000,res=300)
boxplot(diversity_long_gDNA$clones_gDNA~diversity_long_gDNA$recipient.Race,col=COLOR,main="clones by Donor.Source")
dev.off()
summary(glm(clones_gDNA~Donor.Source,data=diversity_long_gDNA))
summary(glm(clones_gDNA~clin,data=diversity_long_gDNA))

###Association with immunosupression
tiff("plot_violin_clones_gDNA_Donor_Source.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==0),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time==0)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==0)],fill=Donor.Source)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=COLOR) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 0",x="Donor.Source", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==6),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time==6)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==6)],fill=Donor.Source)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=COLOR) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 6",x="Donor.Source", y = "Clones")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==24),], 
            aes(factor(diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time==24)]), 
                diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==24)],fill=Donor.Source)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=COLOR) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") + labs(title="time 24",x="Donor.Source", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==6)] ~ diversity_long_gDNA$Donor.Source[which(diversity_long_gDNA$time==6)]))

####Longitudinal model for Donor.Source by time
fm_null <- lmer(clones_est_gDNA ~  time +Donor.Source + (time | Sample_id),
                data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_est_gDNA ~  time*Donor.Source + (time | Sample_id) ,
                data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_donor_source_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_gDNA, colour = Donor.Source)) +
  scale_colour_manual(values=COLOR) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int2_clones_donor_source_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("Donor.Source", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=Donor.Source)) +
  geom_line(aes(time,fit,group=Donor.Source,colour=Donor.Source))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=COLOR)
dev.off()


###Adjust the association with clin by Donor.Source
fm_null <- lmer(clones_est_gDNA ~ clin + time + Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_est_gDNA ~  clin*time + Donor.Source + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_clin_by_donor_source_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int2_clones_clin_by_donor_source_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()


########
## 3. Recipient Age 
########
tiff("boxplot_Rec_age_clones_gDNA.tif",h=2000,w=2000,res=300)
plot(diversity_long_gDNA$clones_gDNA~diversity_long_gDNA$Recipient.Age.when.had.Tx,main="clones by Recipient age")
dev.off()
summary(glm(clones_gDNA~Recipient.Age.when.had.Tx,data=diversity_long_gDNA))
summary(glm(clones_gDNA~clin,data=diversity_long_gDNA))

###Association with RecAge
tiff("plot_violin_clones_gDNA_Rec_Age.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==0),],aes(diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time==0)], 
            diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==0)]))+geom_point() + geom_smooth(method='lm') + labs(title="time 0",x="Rec.Age", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==6),],aes(diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time==6)], 
           diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==6)]))+geom_point() + geom_smooth(method='lm') + labs(title="time 6",x="Rec.Age", y = "Clones")
p3 =ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time==24),],aes(diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time==24)], 
           diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==24)]))+geom_point() + geom_smooth(method='lm') + labs(title="time 24",x="Rec.Age", y = "Clones")

multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time==24)] ~ diversity_long_gDNA$Recipient.Age.when.had.Tx[which(diversity_long_gDNA$time==24)]))

####Longitudinal model for RecAge by time
fm_null <- lmer(clones_est_gDNA ~  time +Recipient.Age.when.had.Tx + (time | Sample_id),
                data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_est_gDNA ~  time*Recipient.Age.when.had.Tx + (time | Sample_id) ,
                data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 


###Adjust the association with clin by Donor.Source
fm_null <- lmer(clones_est_gDNA ~ clin + time + Recipient.Age.when.had.Tx + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(clones_est_gDNA ~  clin*time + Recipient.Age.when.had.Tx + (time | Sample_id),data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_clin_by_recAge_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int2_clones_clin_by_recAge_gDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,group=clin,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()




###########
### AR  ###
###########
diversity_AR<-diversity[which(diversity$clin=="AR" | diversity$clin=="pre-AR"),]
diversity_AR_gDNA<-diversity_AR[which(diversity_AR$reads_gDNA>=100),]
diversity_AR_gDNA$clin<-relevel(diversity_AR_gDNA$clin,ref="pre-AR")

##Number of clones
ggplot(diversity_AR_gDNA, aes(clin, clones_gDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Clones")

##Entropy
ggplot(diversity_AR_gDNA, aes(clin, entropy_gDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Entropy")

##################################
#####Analysis by time and clin ##
##################################
tiff("plot_violin_entropy_est_AR_gDNA.tiff",h=2000,w=1800,res=300)
ggplot(diversity_AR_gDNA, aes(clin,entropy_est_gDNA,fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("goldenrod","firebrick3")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange") + labs(x="Clin", y = "Clones")
dev.off()
summary(glm(diversity_AR_gDNA$entropy_est_gDNA ~ diversity_AR_gDNA$clin))


##Fitting a longitudinal model 
##Clones
fm_null <- lmer(diversity_AR_gDNA$entropy_est_gDNA~ 1 + (1 | Sample_id),data=diversity_AR_gDNA)
fm_full <- lmer(diversity_AR_gDNA$entropy_est_gDNA ~  clin + (1 | Sample_id) ,data=diversity_AR_gDNA)
anova(fm_null,fm_full)

tiff("plot_lmer_entropy_AR_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = clin, y = diversity_AR_gDNA$entropy_est_gDNA, colour = clin)) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  geom_point(size=3) +
  geom_line(aes(group = Sample_id),size=0.7,color="black") 
print(p)
dev.off()

tiff("Plot_lmer2_entropy_AR_gDNA.tiff",h=1700,w=2000,res=300)
fill<-c("goldenrod","firebrick3")
plot(Effect(c("clin"), fm_full))
dev.off()
















#############
### cDNA ####
#############
cor(diversity$clones_cDNA,diversity$clones_est_cDNA,use = "complete.obs")
tiff("plot_cor_clones_cDNA.tiff",res=300,w=1500,h=1500)
qplot(clones_cDNA, clones_est_cDNA, data = diversity, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity$entropy_cDNA,diversity$entropy_est_cDNA,use = "complete.obs")
tiff("plot_cor_entropy_cDNA.tiff",res=300,w=1500,h=1500)
qplot(entropy_cDNA, entropy_est_cDNA, data = diversity, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()

cor(diversity$simpson_cDNA,diversity$simpson_est_cDNA,use = "complete.obs")
tiff("plot_cor_simpson_cDNA.tiff",res=300,w=1500,h=1500)
qplot(simpson_cDNA, simpson_est_cDNA, data = diversity, colour = clin) +
  scale_colour_manual(values=c("firebrick3","chartreuse4","dodgerblue3","darkorange2","goldenrod"))
dev.off()


###Longitudinal Data
diversity_long<-diversity[which(diversity$clin=="NP" | diversity$clin=="PNR" | diversity$clin=="PR"),]
diversity_long_cDNA<-diversity_long[which(diversity_long$reads_cDNA>=100),]
table(diversity_long_cDNA$clin,diversity_long_cDNA$time)

p1<-ggplot(data=diversity_long_cDNA, aes(x=time, y=clones_est_cDNA, group=Sample_id, shape=clinLong, color=clinLong)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal cDNA")

##Number of clones
g1<-ggplot(diversity_long_cDNA[1:24,], aes(time, clones_cDNA)) + scale_y_continuous(limit = c(0,5000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_cDNA[25:48,], aes(time, clones_cDNA)) + geom_point() + scale_y_continuous(limit = c(0,5000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_cDNA[49:69,], aes(time, clones_cDNA)) + geom_point() + scale_y_continuous(limit = c(0,5000)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
multiplot(g1, g2, g3, rows=3)

##Entropy
g1<-ggplot(diversity_long_cDNA[1:24,], aes(time, entropy_cDNA)) + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_cDNA[25:48,], aes(time, entropy_cDNA)) + geom_point() + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_cDNA[49:69,], aes(time, entropy_cDNA)) + geom_point() + scale_y_continuous(limit = c(6,12)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
multiplot(g1, g2, g3, rows=3)

##Simpson
g1<-ggplot(diversity_long_cDNA[1:24,], aes(time, simpson_cDNA)) + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_point() + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g2<-ggplot(diversity_long_cDNA[25:48,], aes(time, simpson_cDNA)) + geom_point() + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
g3<-ggplot(diversity_long_cDNA[49:69,], aes(time, simpson_cDNA)) + geom_point() + scale_y_continuous(limit = c(0,0.03)) + 
  scale_x_continuous(breaks = c(0,6,12,24,32)) + facet_grid(clin ~ Sample_id) + 
  geom_smooth(method = "lm",se = F, colour = "steelblue", size = 1)+ labs(x = "time", y = "Clones")
multiplot(g1, g2, g3, rows=3)

##################################
#####Analysis by time and clin ##
##################################
tiff("plot_violin_clones_cDNA.tiff",h=2000,w=1800,res=300)
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time==6),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time==6)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time==6)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl,geom="pointrange") + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time==24),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time==24)]), 
                diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time==24)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange") + labs(title="time 24",x="Clin", y = "Clones")

multiplot(p2,p3)
dev.off()
summary(glm(diversity_long_cDNA$clones_cDNA[which(diversity_long_cDNA$time==6)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time==6)]))

tiff("plot_violin_entropy_cDNA.tiff",h=2000,w=1800,res=300)
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time==6),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time==6)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time==6)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange") + labs(title="time 6",x="Clin", y = "Entropy")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time==24),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time==24)]), 
                diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time==24)],fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange") + labs(title="time 24",x="Clin", y = "Entropy")

multiplot(p2,p3)
dev.off()

summary(glm(diversity_long_cDNA$entropy_cDNA[which(diversity_long_cDNA$time==6)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time==6)]))


########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
fm_null <- lmer(diversity_long_cDNA$clones_est_cDNA ~ clin + time + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(diversity_long_cDNA$clones_est_cDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

#ggplot(fm_full, aes(time, diversity_long_cDNA$clones_cDNA, color=clin)) +
#  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
#  stat_summary(fun.data=mean_se, geom="pointrange") +
# stat_summary(aes(y=.fitted), fun.y=mean, geom="line")
#  geom_line(aes(y = predict(fm_full))) 

tiff("plot_lmer_clones_est_cDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_cDNA$clones_est_cDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int1_clones_est_cDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_clones_est_cDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper,group=clin),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()

##Diversity
fm_null <- lmer(diversity_long_cDNA$entropy_est_cDNA ~ clin + time + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(diversity_long_cDNA$entropy_est_cDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_entropy_est_cDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_cDNA$entropy_est_cDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_line(aes(y = predict(fm_full)),size=1) 
print(p)
dev.off()

tiff("Plot_int1_entropy_est_cDNA.tiff",h=1000,w=2000,res=300)
fill<-c("chartreuse4", "dodgerblue3","darkorange2")
plot(Effect(c("clin", "time"), fm_full))
dev.off()

tiff("Plot_int2_entropy_est_cDNA.tiff",h=1000,w=1500,res=300)
eff_df <- data.frame(Effect(c("clin", "time"), fm_full))
ggplot(eff_df,aes(time,fit,group=clin)) +
  geom_line(aes(time,fit,colour=clin))+
  geom_ribbon(aes(time,ymin=lower,ymax=upper,group=clin),alpha=0.3) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
dev.off()


###########
### AR  ###
###########
diversity_AR<-diversity[which(diversity$clin=="AR" | diversity$clin=="pre-AR"),]
diversity_AR_cDNA<-diversity_AR[which(diversity_AR$reads_cDNA>=100),]
diversity_AR_cDNA$clin<-relevel(diversity_AR_cDNA$clin,ref="pre-AR")

##Number of clones
ggplot(diversity_AR_cDNA, aes(clin, clones_cDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Clones")

##Entropy
ggplot(diversity_AR_cDNA, aes(clin, entropy_cDNA,group=Sample_id)) +
  geom_point() + geom_line(color="firebrick3") + facet_grid(~Sample_id) + 
  labs(x = "clin", y = "Entropy")

##################################
#####Analysis by time and clin ##
##################################
tiff("plot_violin_clone_est_AR_cDNA.tiff",h=2000,w=1800,res=300)
ggplot(diversity_AR_cDNA, aes(clin,clones_est_cDNA,fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("goldenrod","firebrick3")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange") + labs(x="Clin", y = "Clones")
dev.off()
summary(glm(diversity_AR_cDNA$clones_est_cDNA ~ diversity_AR_cDNA$clin))

tiff("plot_violin_entropy_AR_cDNA.tiff",h=2000,w=1800,res=300)
ggplot(diversity_AR_cDNA, aes(clin,entropy_cDNA,fill=clin)) + 
  geom_violin(trim=FALSE) + scale_fill_manual(values=c("goldenrod","firebrick3")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange") + labs(x="Clin", y = "Entropy")
dev.off()
summary(glm(diversity_AR_cDNA$entropy_est_cDNA ~ diversity_AR_cDNA$clin))


##Fitting a longitudinal model 
##Clones
fm_null <- lmer(diversity_AR_cDNA$entropy_est_cDNA ~ 1 + (1 | Sample_id),data=diversity_AR_cDNA)
fm_full <- lmer(diversity_AR_cDNA$entropy_est_cDNA ~  clin + (1 | Sample_id) ,data=diversity_AR_cDNA)
anova(fm_null,fm_full)

tiff("plot_lmer_entropy_est_AR_cDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = clin, y = diversity_AR_cDNA$entropy_est_cDNA, colour = clin)) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  geom_point(size=3) +
  geom_line(aes(group = Sample_id),size=0.7,color="black") 
print(p)
dev.off()

tiff("Plot_lmer2_entropy_est_AR_cDNA.tiff",h=1700,w=2000,res=300)
fill<-c("goldenrod","firebrick3")
plot(Effect(c("clin"), fm_full))
dev.off()



