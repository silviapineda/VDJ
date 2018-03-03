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
### Date: June, 2017
############################################################################################
library(circlize)
library("RColorBrewer")
library(gtools)
library(lme4)
library("randomForest")
library("VSURF")
library(pheatmap)
library(ggplot2)
library(reshape2)

setwd("/Users/Pinedasans/VDJ/ResultsAllClones/")
load("/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")

############
### cDNA ###
############
data_cDNA<-data_merge[which(data_merge$amplification_template=="cDNA"),]

### longitudinal
data_cDNA_long<-data_cDNA[which(data_cDNA$clin!="AR" & data_cDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]

##Filter for IGHD isotype
data_cDNA_long_IGHD<-data_cDNA_long[which(data_cDNA_long$isotype=="IGHA"),]
#id_clones<-match(data_gDNA_long$V_J_lenghCDR3_CloneId,data_cDNA_long_IGHD$V_J_lenghCDR3_CloneId)
#data_cDNA_long_IGHD_common<-data_cDNA_long_IGHD[na.omit(id_clone),]

reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_IGHA>100),]

id<-match(data_cDNA_long_IGHD$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_cDNA_long_IGHD_qc<-data_cDNA_long_IGHD[which(is.na(id)==F),]
data_cDNA_long_IGHD_qc$specimen_label<-factor(data_cDNA_long_IGHD_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_cDNA_long_IGHD_qc$specimen_label,data_cDNA_long_IGHD_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_qc$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time[id.spec],
              reads_clones_annot_Long_qc$Sample_id[id.spec],reads_clones_annot_Long_qc$Individual.id[id.spec])
colnames(vgenes)[60:63]<-c("clin","time","Sample_id","Individual_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-4))
for (i in 1:59){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_qc$clones_IGHA[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:59]
rownames(vusage)<-vgenes$subject_id

##Delete sample 27 time 12
vgenes2<-vgenes[which(rownames(vgenes)!="M154-S001"),]
vusage2<-vusage[which(rownames(vgenes)!="M154-S001"),]

#Replace some time points
vgenes2$time<-replace(vgenes2$time,vgenes2$time==13,24)
vgenes2$time<-replace(vgenes2$time,vgenes2$time==12,6)

vgenes_filter<-vgenes2
vusage_filter<-vusage2
rownames(vusage_filter)<-rownames(vgenes_filter)

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.05,0)
###Those who are in lesss than 10%
vusage_filter<-vusage_filter[,which(apply(xx,2,function(x) sum(x==0))<=47)]
##39 genes in total
vgenes_filter$clin<-factor(vgenes_filter$clin)

#########
## Statistical analysis to find significant vgenes
########

##ANOVA by time points
p.PNR<-matrix(NA,3,39)
p.PR<-matrix(NA,3,39)
for(i in 1:39){
  #time6
  fit6<-summary(glm(vusage_filter[which(vgenes_filter$time==6),i] ~ vgenes_filter$clin[which(vgenes_filter$time==6)]))
  p.PNR[1,i]<-fit6$coefficients[2,4]
  p.PR[1,i]<-fit6$coefficients[3,4]
  #time24
  fit24<-summary(glm(vusage_filter[which(vgenes_filter$time==24),i] ~ vgenes_filter$clin[which(vgenes_filter$time==24)]))
  p.PNR[2,i]<-fit24$coefficients[2,4]
  p.PR[2,i]<-fit24$coefficients[3,4]
  #all
  fit<-summary(glm(vusage_filter[,i] ~ vgenes_filter$clin))
  p.PNR[3,i]<-fit$coefficients[2,4]
  p.PR[3,i]<-fit$coefficients[3,4]
}
##PNR not statistically significance
colnames(p.PNR)<-colnames(vusage_filter)
colnames(p.PR)<-colnames(vusage_filter)
rownames(p.PNR)<-c("time6","time24","all")
rownames(p.PR)<-c("time6","time24","all")
p.PR<0.05

vgenes_sign_6<-names(which(p.PR[1,]<0.1))
vgenes_sign_24<-names(which(p.PR[2,]<0.1))

names(which(p.adjust(p.PR[1,],method = "fdr")<0.1)) 
names(which(p.adjust(p.PR[2,],method = "fdr")<0.1)) 

id_sign_6<-match(vgenes_sign_6,colnames(vusage_filter))
vusage_sign_6<-vusage_filter[,id_sign_6]
id_sign_24<-match(vgenes_sign_24,colnames(vusage_filter))
vusage_sign_24<-vusage_filter[,id_sign_24]

vusage_sign_6<-vusage_sign_6[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR"),]
vusage_sign_24<-vusage_sign_24[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR"),]

annotation_col_6 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]))
annotation_col_24 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]))


rownames(annotation_col_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]
rownames(annotation_col_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]

rownames(vusage_sign_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]
rownames(vusage_sign_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]

tiff("boxplot_IGHV3-23_IGHG.tiff",res=300,w=1000,h=1500)
#par(mfrow=c(1,3))
#boxplot(vusage_sign_6[,"IGHV3-23"]~annotation_col_6$clin,
#        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 6")
boxplot(vusage_sign_24[,"IGHV3-23"]~annotation_col_24$clin,
        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.8),ylab="IHGV3-23 expression (IGHG)",main="time 24")
dev.off()

