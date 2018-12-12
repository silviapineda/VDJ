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
### gDNA ###
############
data_gDNA<-data_merge[which(data_merge$amplification_template=="gDNA"),]
### longitudinal
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]

#For this analysis we are putting a cut-off on clones because v-genes can be biased at low clonality 
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_gDNA>130),]

id<-match(data_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_gDNA_long_qc$specimen_label,data_gDNA_long_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_qc$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time[id.spec],
              reads_clones_annot_Long_qc$Sample_id[id.spec],reads_clones_annot_Long_qc$Individual.id[id.spec],reads_clones_annot_Long_qc$Cause.of.ESRD[id.spec],
              reads_clones_annot_Long_qc$ESRD[id.spec],reads_clones_annot_Long_qc$immunosuppression[id.spec],
              reads_clones_annot_Long_qc$Donor.Source[id.spec],reads_clones_annot_Long_qc$Recipient.Gender[id.spec],
              reads_clones_annot_Long_qc$Donor.Gender[id.spec],reads_clones_annot_Long_qc$Recipient.Age.when.had.Tx[id.spec],
              reads_clones_annot_Long_qc$Donor.Age[id.spec],reads_clones_annot_Long_qc$recipient.Race[id.spec],
              reads_clones_annot_Long_qc$hla_mismatch[id.spec])
colnames(vgenes)[65:78]<-c("clin","time","Sample_id","Individual_id","causeESRD","ESRD","immunosuppression","Source","RecGender","DonGen","RecAge",
                           "DonAge","RecRace","HLA")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-14))
for (i in 1:64){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:64]
rownames(vusage)<-vgenes$subject_id

##Delete sample 8
vgenes2<-vgenes[which(vgenes$Individual_id!="Individual8"),]
vusage2<-vusage[which(vgenes$Individual_id!="Individual8"),]

#Replace some time points
vgenes2$time<-replace(vgenes2$time,vgenes2$time==13,24)
vgenes2$time<-replace(vgenes2$time,vgenes2$time==12,6)
#vgenes2$time<-replace(vgenes2$time,vgenes2$time==32,24)

vgenes_filter<-vgenes2
vusage_filter<-vusage2
rownames(vusage_filter)<-rownames(vgenes_filter)

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.05,0)
###Those who are in lesss than 10%
vusage_filter<-vusage_filter[,which(apply(xx,2,function(x) sum(x==0))<=57)]
##27 genes in total
vgenes_filter$clin<-factor(vgenes_filter$clin)

NP_usage<-colSums(vusage_filter[which(vgenes_filter$clin=="NP"),])
PNR_usage<-colSums(vusage_filter[which(vgenes_filter$clin=="PNR"),])
PR_usage<-colSums(vusage_filter[which(vgenes_filter$clin=="PR"),])

PR_usage[order(PR_usage)]
PNR_usage[order(PNR_usage)]
NP_usage[order(NP_usage)]

#########
## Statistical analysis to find significant vgenes
########

##ANOVA by time points
p.PNR<-matrix(NA,4,27)
p.PR<-matrix(NA,4,27)
for(i in 1:27){
  #time0
  fit0<-summary(glm(vusage_filter[which(vgenes_filter$time==0),i] ~ vgenes_filter$clin[which(vgenes_filter$time==0)]))
  p.PNR[1,i]<-fit0$coefficients[2,4]
  p.PR[1,i]<-fit0$coefficients[3,4]
  #time6
  fit6<-summary(glm(vusage_filter[which(vgenes_filter$time==6),i] ~ vgenes_filter$clin[which(vgenes_filter$time==6)]))
  p.PNR[2,i]<-fit6$coefficients[2,4]
  p.PR[2,i]<-fit6$coefficients[3,4]
  #time24
  fit24<-summary(glm(vusage_filter[which(vgenes_filter$time==24),i] ~ vgenes_filter$clin[which(vgenes_filter$time==24)]))
  p.PNR[3,i]<-fit24$coefficients[2,4]
  p.PR[3,i]<-fit24$coefficients[3,4]
  #all
  fit<-summary(glm(vusage_filter[,i] ~ vgenes_filter$clin))
  p.PNR[4,i]<-fit$coefficients[2,4]
  p.PR[4,i]<-fit$coefficients[3,4]
}
##PNR not statistically significance
colnames(p.PNR)<-colnames(vusage_filter)
colnames(p.PR)<-colnames(vusage_filter)
rownames(p.PNR)<-c("time0","time6","time24","all")
rownames(p.PR)<-c("time0","time6","time24","all")
p.PR<0.05

vgenes_sign_0<-names(which(p.PR[1,]<0.05))
vgenes_sign_6<-names(which(p.PR[2,]<0.05))
vgenes_sign_24<-names(which(p.PR[3,]<0.05))

names(which(p.adjust(p.PR[1,],method = "fdr")<0.1)) 
names(which(p.adjust(p.PR[2,],method = "fdr")<0.1)) 
names(which(p.adjust(p.PR[3,],method = "fdr")<0.1)) 

##Heatmap with significant results
id_sign_0<-match(vgenes_sign_0,colnames(vusage_filter))
vusage_sign_0<-vusage_filter[,id_sign_0]
id_sign_6<-match(vgenes_sign_6,colnames(vusage_filter))
vusage_sign_6<-vusage_filter[,id_sign_6]
id_sign_24<-match(vgenes_sign_24,colnames(vusage_filter))
vusage_sign_24<-vusage_filter[,id_sign_24]

vusage_sign_0<-vusage_sign_0[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR"),]
vusage_sign_6<-vusage_sign_6[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR"),]
vusage_sign_24<-vusage_sign_24[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR"),]

annotation_col_0 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
  causeESRD = factor(vgenes_filter$causeESRD[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
  immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
  source = factor(vgenes_filter$Source[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]))
annotation_col_6 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
  causeESRD = factor(vgenes_filter$causeESRD[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
  immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
  source = factor(vgenes_filter$Source[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]))
annotation_col_24 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
  causeESRD = factor(vgenes_filter$causeESRD[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
  immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
  source = factor(vgenes_filter$Source[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]))


rownames(annotation_col_0)<-vgenes_filter$Individual_id[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]
rownames(annotation_col_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]
rownames(annotation_col_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]

rownames(vusage_sign_0)<-vgenes_filter$Individual_id[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]
rownames(vusage_sign_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]
rownames(vusage_sign_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]



fill=c("chartreuse4","darkorange2")
fill2=brewer.pal(7,"Set2")
fill3 = brewer.pal(3,"Set1")
fill4 = b=brewer.pal(3,"Set3")

ann_colors = list (clin = c("NP" = fill[1], "PR" = fill[2]),
                   causeESRD = c("Cortical Necrosis" = fill2[1], "Focal Segmental Glomerulosclerosis" = fill2[2], 
                                 "Obstructive Uropathy" = fill2[3], "Other" =fill2[4], "Reflux Nephropathy" = fill2[5],
                                 "Unknown" =fill2[6], "Polycystic Kidney Disease" =fill2[7]))
                   #immunosuppression = c("Steroid-free" = fill3[1], "Steroid-based" = fill3[2]),
                   #source = c("Cadaver" = fill4[1],"Living/related"=fill4[2]))

colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_0.tiff",res=300,w=2100,h=1200)
pheatmap(t(vusage_sign_0),annotation_col = annotation_col_0[,1:2],fontsize = 8,main="time 0"
         ,annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()
tiff("heatmap_vgene_6.tiff",res=300,w=2100,h=1200)
pheatmap(t(vusage_sign_6),annotation_col = annotation_col_6 [1:2],fontsize = 8,main="time 6",
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()
tiff("heatmap_vgene_24.tiff",res=300,w=2100,h=1200)
pheatmap(t(vusage_sign_24),annotation_col = annotation_col_24[1:2], fontsize = 8, main="time 24",
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

tiff("boxplot_IGHV3-23.tiff",res=300,w=1500,h=1000)
par(mfrow=c(1,3))
boxplot(vusage_sign_0[,"IGHV3-23"]~annotation_col_0$clin,
        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 0")
boxplot(vusage_sign_6[,"IGHV3-23"]~annotation_col_6$clin,
        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 6")
boxplot(vusage_sign_24[,"IGHV3-23"]~annotation_col_24$clin,
        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 24")
dev.off()

tiff("boxplot_IGHV3-23_ESRD.tiff",res=300,w=2500,h=2000)
par(mfrow=c(1,3))
boxplot(vusage_sign_0[,"IGHV3-23"]~annotation_col_0$ESRD,
        ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 0",las=2,cex.axis=0.8)
boxplot(vusage_sign_6[,"IGHV3-23"]~annotation_col_6$ESRD,
        ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 6",las=2,cex.axis=0.8)
boxplot(vusage_sign_24[,"IGHV3-23"]~annotation_col_24$ESRD,
       ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 24",las=2,cex.axis=0.8)
dev.off()

summary(glm(vusage_sign_0[,"IGHV3-23"]~annotation_col_0$clin + annotation_col_0$ESRD))
summary(glm(vusage_sign_6[,"IGHV3-23"]~annotation_col_6$clin + annotation_col_6$ESRD))
summary(glm(vusage_sign_24[,"IGHV3-23"]~annotation_col_24$clin + annotation_col_24$ESRD))

summary(glm(vusage_sign_0[,"IGHV3-23"]~factor(annotation_col_0$ESRD)))
summary(glm(vusage_sign_6[,"IGHV3-23"]~factor(annotation_col_6$ESRD)))
summary(glm(vusage_sign_24[,"IGHV3-23"]~factor(annotation_col_24$ESRD)))

COLOR=c("chartreuse4","darkorange2")
tiff("boxplot_clones_causeESRD_time0.tiff",h=2000,w=3500,res=300)
ggplot(annotation_col_0,aes(factor(causeESRD),vusage_sign_0[,"IGHV3-23"],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 0",x="ESRD", y = "IGHV3-23 expression")
dev.off()
tiff("boxplot_clones_causeESRD_time6.tiff",h=2000,w=3500,res=300)
ggplot(annotation_col_6,aes(factor(causeESRD),vusage_sign_6[,"IGHV3-23"],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 6",x="ESRD", y = "IGHV3-23 expression")
dev.off()
tiff("boxplot_clones_causeESRD_time24.tiff",h=2000,w=3500,res=300)
ggplot(annotation_col_24,aes(factor(causeESRD),vusage_sign_24[,"IGHV3-23"],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time 24",x="ESRD", y = "IGHV3-23 expression")
dev.off()


#####Considering the 3 clinical outcomes
vusage_sign_0<-vusage_sign_0[which(vgenes_filter$time==0),]
vusage_sign_6<-vusage_sign_6[which(vgenes_filter$time==6),]
vusage_sign_24<-vusage_sign_24[which(vgenes_filter$time==24),]

annotation_col_0 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==0)]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==0)]))
#immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
#source = factor(vgenes_filter$Source[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]))
annotation_col_6 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==6)]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==6)]))
#immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
#source = factor(vgenes_filter$Source[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]))
annotation_col_24 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==24)]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==24)]))
#immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
#source = factor(vgenes_filter$Source[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]))


rownames(annotation_col_0)<-vgenes_filter$Individual_id[which(vgenes_filter$time==0)]
rownames(annotation_col_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6)]
rownames(annotation_col_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24)]

rownames(vusage_sign_0)<-vgenes_filter$Individual_id[which(vgenes_filter$time==0)]
rownames(vusage_sign_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6)]
rownames(vusage_sign_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24)]



fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set2")
#fill3 = brewer.pal(3,"Set1")
#fill4 = b=brewer.pal(3,"Set3")

ann_colors = list (clin = c("NP" = fill[1],"PNR"= fill[2],"PR" = fill[3]),
                   ESRD = c("Other/Uk" = fill2[1], "Reflux" = fill2[2], "Non-Immune/Strutural" = fill2[3]))
#immunosuppression = c("Steroid-free" = fill3[1], "Steroid-based" = fill3[2]),
#source = c("Cadaver" = fill4[1],"Living/related"=fill4[2]))

colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_0.tiff",res=300,w=1500,h=1200)
pheatmap(t(vusage_sign_0),annotation_col = annotation_col_0,fontsize = 8,main="time 0"
         ,annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()
tiff("heatmap_vgene_6.tiff",res=300,w=1500,h=1200)
pheatmap(t(vusage_sign_6),annotation_col = annotation_col_6,fontsize = 8,main="time 6",
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()
tiff("heatmap_vgene_24.tiff",res=300,w=1500,h=1200)
pheatmap(t(vusage_sign_24),annotation_col = annotation_col_24, fontsize = 8, main="time 24",
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

tiff("boxplot_IGHV3-23.tiff",res=300,w=1500,h=1000)
par(mfrow=c(1,3))
boxplot(vusage_sign_0[,"IGHV3-23"]~annotation_col_0$clin,
        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 0")
boxplot(vusage_sign_6[,"IGHV3-23"]~annotation_col_6$clin,
        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 6")
boxplot(vusage_sign_24[,"IGHV3-23"]~annotation_col_24$clin,
        col=c("chartreuse4","darkorange2"),ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 24")
dev.off()

tiff("boxplot_IGHV3-23_ESRD.tiff",res=300,w=1500,h=1000)
par(mfrow=c(1,3))
boxplot(vusage_sign_0[,"IGHV3-23"]~annotation_col_0$ESRD,
        ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 0")
boxplot(vusage_sign_6[,"IGHV3-23"]~annotation_col_6$ESRD,
        ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 6")
boxplot(vusage_sign_24[,"IGHV3-23"]~annotation_col_24$ESRD,
        ylim=c(0.0,0.5),ylab="IHGV3-23 expression",main="time 24")
dev.off()


summary(glm(vusage_filter[which(vgenes_filter$time==0),"IGHV3-23"]~vgenes_filter$clin[which(vgenes_filter$time==0)]))
summary(glm(vusage_filter[which(vgenes_filter$time==6),"IGHV3-23"]~vgenes_filter$clin[which(vgenes_filter$time==6)]))
summary(glm(vusage_filter[which(vgenes_filter$time==24),"IGHV3-23"]~vgenes_filter$clin[which(vgenes_filter$time==24)]))

summary(glm(vusage_filter[which(vgenes_filter$time==0),"IGHV3-23"]~vgenes_filter$clin[which(vgenes_filter$time==0)]
            +vgenes_filter$HLA[which(vgenes_filter$time==0)]))

summary(glm(vusage_filter[which(vgenes_filter$time==6),"IGHV3-23"]~vgenes_filter$clin[which(vgenes_filter$time==6)]
            +vgenes_filter$HLA[which(vgenes_filter$time==6)]))

summary(glm(vusage_filter[which(vgenes_filter$time==24),"IGHV3-23"]~vgenes_filter$clin[which(vgenes_filter$time==24)]
              +vgenes_filter$HLA[which(vgenes_filter$time==24)]))

summary(glm(vusage_filter[which(vgenes_filter$time==0),"IGHV3-23"]~vgenes_filter$ESRD[which(vgenes_filter$time==0)]))
summary(glm(vusage_filter[which(vgenes_filter$time==6),"IGHV3-23"]~vgenes_filter$ESRD[which(vgenes_filter$time==6)]))
summary(glm(vusage_filter[which(vgenes_filter$time==24),"IGHV3-23"]~vgenes_filter$ESRD[which(vgenes_filter$time==24)]))


###Longitudinal analysis
p<-NULL
for(i in 1:27){
  fm_null <- lmer(vusage_filter[,i] ~ vgenes_filter$clin + vgenes_filter$time + (vgenes_filter$time | vgenes_filter$Sample_id),REML = F)
  fm_full <- lmer(vusage_filter[,i] ~  vgenes_filter$clin*vgenes_filter$time + (vgenes_filter$time | vgenes_filter$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##
} 
names(p)<-colnames(vusage_filter)
sign_id<-which(p<0.05)
for(i in sign_id){
  fm_null <- lmer(vusage_filter[,i] ~ clin + time + (time | Individual_id),data=vgenes_filter,REML = F)
  fm_full <- lmer(vusage_filter[,i] ~ clin*time + (time | Individual_id),data=vgenes_filter,REML = F)
  tiff(paste0("plot_lmer_gDNA",colnames(vusage_filter)[i],".tiff"),h=1200,w=1400,res=300)
  plot <- ggplot(fm_full, aes(x = time , y = vusage_filter[,i], colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=1.2) +
    geom_smooth(method="lm",size=0.8) +
    labs(x = "Time (months)",y = paste0(colnames(vusage_filter)[i]," expression")) + theme_bw() + theme_light()
  
  print(plot)
  dev.off()
} 

#################################
####### V  J junction #########
###############################
data_gDNA_long$vjgenes<- paste(data_gDNA_long$v_gene,data_gDNA_long$j_gene,sep="_")
vjgenes<-as.data.frame(unclass(table(data_gDNA_long$specimen_label,data_gDNA_long$vjgenes)))
id.spec<-match(rownames(vjgenes),reads_clones_annot_Long_qc$specimen_id)
vjgenes<-cbind(vjgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time[id.spec],
               reads_clones_annot_Long_qc$Sample_id[id.spec],reads_clones_annot_Long_qc$Individual.id[id.spec])
colnames(vjgenes)[327:330]<-c("clin","time","Sample_id","individual_id")

vjusage<-matrix(NA,nrow(vjgenes),(ncol(vjgenes)-4))
for (i in 1:326){
  vjusage[,i]<-vjgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(vjusage)<-colnames(vjgenes)[1:326]
rownames(vjusage)<-vjgenes$subject_id

##Delete sample 8
vjgenes2<-vjgenes[which(vjgenes$individual_id!="Indiviudal8"),]
vjusage2<-vjusage[which(vjgenes$individual_id!="Indiviudal8"),]

#Replace some time points
vjgenes2$time<-replace(vjgenes2$time,vjgenes2$time==13,24)
vjgenes2$time<-replace(vjgenes2$time,vjgenes2$time==12,6)
vjgenes2$time<-replace(vjgenes2$time,vjgenes2$time==32,24)

vjgenes_filter<-vjgenes2
vjusage_filter<-vjusage2
rownames(vjusage_filter)<-rownames(vjgenes_filter)

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vjusage_filter,vjusage_filter<0.01,0)
###Those who are in lesss than 10%
vjusage_filter<-vjusage_filter[,which(apply(xx,2,function(x) sum(x==0))<=57)]
##112 genes in total
vjgenes_filter$clin<-factor(vjgenes_filter$clin)


#########
## Statistical analysis to find significant vgenes
########

##ANOVA by time points
p.PNR<-matrix(NA,4,112)
p.PR<-matrix(NA,4,112)
for(i in 1:112){
  #time0
  fit0<-summary(glm(vjusage_filter[which(vjgenes_filter$time==0),i] ~ vjgenes_filter$clin[which(vjgenes_filter$time==0)]))
  p.PNR[1,i]<-fit0$coefficients[2,4]
  p.PR[1,i]<-fit0$coefficients[3,4]
  #time6
  fit6<-summary(glm(vjusage_filter[which(vjgenes_filter$time==6),i] ~ vjgenes_filter$clin[which(vjgenes_filter$time==6)]))
  p.PNR[2,i]<-fit6$coefficients[2,4]
  p.PR[2,i]<-fit6$coefficients[3,4]
  #time24
  fit24<-summary(glm(vjusage_filter[which(vjgenes_filter$time==24),i] ~ vjgenes_filter$clin[which(vjgenes_filter$time==24)]))
  p.PNR[3,i]<-fit24$coefficients[2,4]
  p.PR[3,i]<-fit24$coefficients[3,4]
  #all
  fit<-summary(glm(vjusage_filter[,i] ~ vjgenes_filter$clin))
  p.PNR[4,i]<-fit$coefficients[2,4]
  p.PR[4,i]<-fit$coefficients[3,4]
}
##PNR not statistically significance
colnames(p.PNR)<-colnames(vjusage_filter)
colnames(p.PR)<-colnames(vjusage_filter)
rownames(p.PNR)<-c("time0","time6","time24","all")
rownames(p.PR)<-c("time0","time6","time24","all")
p.PR<0.05

vjgenes_sign_0<-names(which(p.PR[1,]<0.05))
vjgenes_sign_6<-names(which(p.PR[2,]<0.05))
vjgenes_sign_24<-names(which(p.PR[3,]<0.05))

p.adjust(p.PR[1,],method = "fdr")
p.adjust(p.PR[2,],method = "fdr")
p.adjust(p.PR[3,],method = "fdr")


##Heatmap with significant results
id_sign_0<-match(vjgenes_sign_0,colnames(vjusage_filter))
vjusage_sign_0<-vjusage_filter[,id_sign_0]
id_sign_6<-match(vjgenes_sign_6,colnames(vjusage_filter))
vjusage_sign_6<-vjusage_filter[,id_sign_6]
id_sign_24<-match(vjgenes_sign_24,colnames(vjusage_filter))
vjusage_sign_24<-vjusage_filter[,id_sign_24]

vjusage_sign_0<-vjusage_sign_0[which(vjgenes_filter$time==0 & vjgenes_filter$clin!="PNR"),]
vjusage_sign_6<-vjusage_sign_6[which(vjgenes_filter$time==6 & vjgenes_filter$clin!="PNR"),]
vjusage_sign_24<-vjusage_sign_24[which(vjgenes_filter$time==24 & vjgenes_filter$clin!="PNR"),]

annotation_col_0 = data.frame(
  clin = factor(vjgenes_filter$clin[which(vjgenes_filter$time==0 & vjgenes_filter$clin!="PNR")]))
annotation_col_6 = data.frame(
  clin = factor(vjgenes_filter$clin[which(vjgenes_filter$time==6 & vjgenes_filter$clin!="PNR")]))
annotation_col_24 = data.frame(
  clin = factor(vjgenes_filter$clin[which(vjgenes_filter$time==24 & vjgenes_filter$clin!="PNR")]))


rownames(annotation_col_0)<-rownames(vjusage_sign_0)
rownames(annotation_col_6)<-rownames(vjusage_sign_6)
rownames(annotation_col_24)<-rownames(vjusage_sign_24)

fill=c("chartreuse4","darkorange2")
ann_colors = list (clin = c("NP" = fill[1], "PR" = fill[2]))

colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vjgene_0.tiff",res=300,w=3500,h=2500)
pheatmap(t(vjusage_sign_0),cex=1.0,annotation_col = annotation_col_0,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()
tiff("heatmap_vjgene_6.tiff",res=300,w=3500,h=2500)
pheatmap(t(vjusage_sign_6),cex=1.0,annotation_col = annotation_col_6,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()
tiff("heatmap_vjgene_24.tiff",res=300,w=3500,h=2500)
pheatmap(t(vjusage_sign_24),cex=1.0,annotation_col = annotation_col_24, 
         clustering_distance_cols="euclidean",
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()


###Longitudinal analysis
p<-NULL
for(i in 1:112){
  print(i)
  fm_null <- lmer(vjusage_filter[,i] ~ vjgenes_filter$clin + vjgenes_filter$time + (vjgenes_filter$time | vjgenes_filter$Sample_id),REML = F)
  fm_full <- lmer(vjusage_filter[,i] ~  vjgenes_filter$clin*vjgenes_filter$time + (vjgenes_filter$time | vjgenes_filter$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##
} 
names(p)<-colnames(vjusage_filter)
sign_id<-which(p<0.05)
for(i in sign_id){
  fm_null <- lmer(vjusage_filter[,i] ~ clin + time + (time | individual_id),data=vjgenes_filter,REML = F)
  fm_full <- lmer(vjusage_filter[,i] ~ clin*time + (time | individual_id),data=vjgenes_filter,REML = F)
  tiff(paste0("plot_lmer_gDNA_vj",colnames(vjusage_filter)[i],".tiff"),h=1200,w=1400,res=300)
  plot <- ggplot(fm_full, aes(x = time , y = vjusage_filter[,i], colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=1.2) +
    geom_smooth(method="lm",size=0.8) +
    labs(x = "Time (months)",y = colnames(vjusage_filter)[i]) + theme_bw() + theme_light()
  
  print(plot)
  dev.off()
} 



###############
####One circos plot with the median values for each clinical outcome
##############
vjusage_sign_NP<-vjusage_sign[which(vjgenes_filter$clin=="NP"),]
vjusage_sign_PNR<-vjusage_sign[which(vjgenes_filter$clin=="PNR"),]
vjusage_sign_PR<-vjusage_sign[which(vjgenes_filter$clin=="PR"),]
vjusage_sign_NP_mean<-apply(vjusage_sign_NP,2,mean)
vjusage_sign_PNR_mean<-apply(vjusage_sign_PNR,2,mean)
vjusage_sign_PR_mean<-apply(vjusage_sign_PR,2,mean)

##NP 
splitop<-strsplit(names(vjusage_sign_NP_mean),"_")
v<-unlist(lapply(splitop, `[[`, 1))
j<-unlist(lapply(splitop, `[[`, 2))
vj_usage<-data.frame(cbind(v,j,as.numeric(vjusage_sign_NP_mean)))
colnames(vj_usage)<-c("link1","link2","count")
vj_usage_df<-vj_usage
vj_usage_df$count<-as.numeric(as.character(vj_usage_df$count))

circos.clear()
tiff(paste("circos_VJ_NP.tiff",sep=""),res=300,h=1800,w=1800)
circos.par(gap.degree = c(rep(2, length(unique(v))-1), 10,
                          rep(2, length(unique(j))-1), 10))

grid.col<-c(colorRampPalette(brewer.pal(6,"Reds"))(length(unique(v))),
            colorRampPalette(brewer.pal(6,"Purples"))(length(unique(j))))
names(grid.col)<-c(as.character(unique(v)),as.character(unique(j)))

chordDiagram(vj_usage_df, grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
}, bg.border = NA) 

dev.off()

##PNR 
splitop<-strsplit(names(vjusage_sign_PNR_mean),"_")
v<-unlist(lapply(splitop, `[[`, 1))
j<-unlist(lapply(splitop, `[[`, 2))
vj_usage<-data.frame(cbind(v,j,as.numeric(vjusage_sign_PNR_mean)))
colnames(vj_usage)<-c("link1","link2","count")
vj_usage_df<-vj_usage
vj_usage_df$count<-as.numeric(as.character(vj_usage_df$count))
circos.clear()
tiff(paste("circos_vj_PNR.tiff",sep=""),res=300,h=1800,w=1800)
circos.par(gap.degree = c(rep(2, length(unique(v))-1), 10,
                          rep(2, length(unique(j))-1), 10))

grid.col<-c(colorRampPalette(brewer.pal(9,"Reds"))(length(unique(v))),
            colorRampPalette(brewer.pal(6,"Purples"))(length(unique(j))))
names(grid.col)<-c(as.character(unique(v)),as.character(unique(j)))

chordDiagram(vj_usage_df, grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
}, bg.border = NA) 

dev.off()

##PR 
splitop<-strsplit(names(vjusage_sign_PR_mean),"_")
v<-unlist(lapply(splitop, `[[`, 1))
j<-unlist(lapply(splitop, `[[`, 2))
vj_usage<-data.frame(cbind(v,j,as.numeric(vjusage_sign_PR_mean)))
colnames(vj_usage)<-c("link1","link2","count")
vj_usage_df<-vj_usage
vj_usage_df$count<-as.numeric(as.character(vj_usage_df$count))
circos.clear()
tiff(paste("circos_vj_PR.tiff",sep=""),res=300,h=1800,w=1800)
circos.par(gap.degree = c(rep(2, length(unique(v))-1), 10,
                          rep(2, length(unique(j))-1), 10))

grid.col<-c(colorRampPalette(brewer.pal(9,"Reds"))(length(unique(v))),
            colorRampPalette(brewer.pal(6,"Purples"))(length(unique(j))))
names(grid.col)<-c(as.character(unique(v)),as.character(unique(j)))

chordDiagram(vj_usage_df, grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
}, bg.border = NA) 

dev.off()



############################################
##Comparing three chor diargams by matrix ##
###########################################


##PR
splitop<-strsplit(names(vjusage_sign_PR_mean),"_")
v<-unlist(lapply(splitop, `[[`, 1))
j<-unlist(lapply(splitop, `[[`, 2))
vj_usage<-data.frame(cbind(v,j,as.numeric(vjusage_sign_PR_mean)))
colnames(vj_usage)<-c("link1","link2","count")
vj_usage_df<-vj_usage
vj_usage_df$count<-as.numeric(as.character(vj_usage_df$count))
vj_usage_mat_PR<-acast(vj_usage_df, link1~link2, value.var="count")
vj_usage_mat_PR2<-replace(vj_usage_mat_PR,is.na(vj_usage_mat_PR)==T,0)
circos.clear()

gap.after = c(rep(2, length(unique(v))-1), 10, rep(2, length(unique(j))-1), 10)

circos.par(gap.degree = c(rep(2, length(unique(v))-1), 10,
                          rep(2, length(unique(j))-1), 10))


circos.par(gap.after = gap.after)

grid.col<-c(colorRampPalette(brewer.pal(9,"Reds"))(length(unique(v))),
            colorRampPalette(brewer.pal(6,"Purples"))(length(unique(j))))
names(grid.col)<-c(rownames(vj_usage_mat_PR2), colnames(vj_usage_mat_PR2))


tiff(paste("v_j_junction_gDNA_long_PR.tiff"),res=300,h=1800,w=1800)

chordDiagram(vj_usage_mat_PR2,  grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
}, bg.border = NA) 
dev.off()

##PNR
splitop<-strsplit(names(vjusage_sign_PNR_mean),"_")
v<-unlist(lapply(splitop, `[[`, 1))
j<-unlist(lapply(splitop, `[[`, 2))
vj_usage<-data.frame(cbind(v,j,as.numeric(vjusage_sign_PNR_mean)))
colnames(vj_usage)<-c("link1","link2","count")
vj_usage_df<-vj_usage
vj_usage_df$count<-as.numeric(as.character(vj_usage_df$count))
vj_usage_mat_PNR<-acast(vj_usage_df, link1~link2, value.var="count")
vj_usage_mat_PNR2<-replace(vj_usage_mat_PNR ,is.na(vj_usage_mat_PNR)==T,0)

circos.clear()

percent = sum(abs(vj_usage_mat_PNR2)) / sum(abs(vj_usage_mat_PR2))
blank.degree = (360 - sum(gap.after)) * (1 - percent)

big.gap = (blank.degree - sum(rep(2, 8)))/2
gap.after = c(rep(2, length(unique(v))-1), big.gap, rep(2, length(unique(j))-1), big.gap)
circos.par(gap.after = gap.after, start.degree = -big.gap/2)

grid.col<-c(colorRampPalette(brewer.pal(9,"Reds"))(length(unique(v))),
            colorRampPalette(brewer.pal(6,"Purples"))(length(unique(j))))

names(grid.col)<-c(rownames(vj_usage_mat_PNR2), colnames(vj_usage_mat_PNR2))

tiff(paste("v_j_junction_gDNA_long_PNR.tiff"),res=300,h=1800,w=1800)

chordDiagram(vj_usage_mat_PNR, grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
}, bg.border = NA) 

dev.off()


##NP
splitop<-strsplit(names(vjusage_sign_NP_mean),"_")
v<-unlist(lapply(splitop, `[[`, 1))
j<-unlist(lapply(splitop, `[[`, 2))
vj_usage<-data.frame(cbind(v,j,as.numeric(vjusage_sign_NP_mean)))
colnames(vj_usage)<-c("link1","link2","count")
vj_usage_df<-vj_usage
vj_usage_df$count<-as.numeric(as.character(vj_usage_df$count))
vj_usage_mat_NP<-acast(vj_usage_df, link1~link2, value.var="count")
vj_usage_mat_NP2<-replace(vj_usage_mat_NP ,is.na(vj_usage_mat_NP)==T,0)

circos.clear()

percent = sum(abs(vj_usage_mat_NP2)) / sum(abs(vj_usage_mat_PR2))
blank.degree = (360 - sum(gap.after)) * (1 - percent)

big.gap = (blank.degree - sum(rep(2, 8)))/2
gap.after = c(rep(2, length(unique(v))-1), big.gap, rep(2, length(unique(j))-1), big.gap)
circos.par(gap.after = gap.after, start.degree = -big.gap/2)

grid.col<-c(colorRampPalette(brewer.pal(9,"Reds"))(length(unique(v))),
            colorRampPalette(brewer.pal(6,"Purples"))(length(unique(j))))
names(grid.col)<-c(rownames(vj_usage_mat_NP2), colnames(vj_usage_mat_NP2))

tiff(paste("v_j_junction_gDNA_long_NP.tiff"),res=300,h=1800,w=1800)

chordDiagram(vj_usage_mat_NP, grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
}, bg.border = NA) 

dev.off()






