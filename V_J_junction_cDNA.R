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

#For this analysis we are putting a cut-off on clones because v-genes can be biased at low clonality 
##########################
###Study by isotypes######
##########################

############
# 1. IGHD
###########
data_cDNA_long_IGHD<-data_cDNA_long[which(data_cDNA_long$isotype=="IGHD"),]
reads_clones_annot_Long_IGHD<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_IGHD>100),]
id<-match(data_cDNA_long_IGHD$specimen_label,reads_clones_annot_Long_IGHD$specimen_id)

data_cDNA_long_IGHD_qc<-data_cDNA_long_IGHD[which(is.na(id)==F),]
data_cDNA_long_IGHD_qc$specimen_label<-factor(data_cDNA_long_IGHD_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_cDNA_long_IGHD_qc$specimen_label,data_cDNA_long_IGHD_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_IGHD$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_IGHD$clin[id.spec],reads_clones_annot_Long_IGHD$time[id.spec],
              reads_clones_annot_Long_IGHD$Sample_id[id.spec],reads_clones_annot_Long_IGHD$subject_id[id.spec])
colnames(vgenes)[63:66]<-c("clin","time","Sample_id","subject_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-4))
for (i in 1:62){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_IGHD$clones_IGHD[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:62]
rownames(vusage)<-vgenes$subject_id


#Filter by time
##First manually convert the time points
vgenes$time<-replace(vgenes$time,vgenes$time==32,24)
vgenes$time<-replace(vgenes$time,vgenes$time==13,24)
vgenes$time<-replace(vgenes$time,vgenes$time==12,6)

##Filter for those that have not time 6 or 24
vgenes_filter<-vgenes[which(vgenes$time==6 | vgenes$time==24),]
vusage_filter<-vusage[which(vgenes$time==6 | vgenes$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.002,0)
###Those genes who are in lesss than 10%
vusage_filter2<-vusage_filter[,which(apply(xx,2,function(x) sum(x==0))<50)]
##47 genes in total

##V-genes found in the gDNA Analysis
# "IGHV4.34" "IGHV4.31" "IGHV3.11" "IGHV3.23" "IGHV3.9"  "IGHV3.33" "IGHV3.21" "IGHV1.46" "IGHV3.7"  "IGHV3.13" "IGHV3.53" "IGHV4.39"
# "IGHV1.24" "IGHV5.a"  "IGHV3.20" "IGHV4.59" "IGHV1.f" 

#########
## Statistical analysis to find significant vgenes
########
set.seed(112233)
vusage_time<-data.frame(cbind(vusage_filter2,vgenes_filter$time))
fit<-VSURF(x = vusage_time, y=factor(vgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vusage_time[,fit$varselect.interp])
set.seed(1000)
rf_output <- randomForest(factor(vgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vgene_cDNA_IGHD_RF.tiff",res=300,w=2000,h=2000)
vgenes_filter$timeplot<-ifelse(vgenes_filter$time==0,1,
                               ifelse(vgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vgenes_filter$clin),cex=vgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "IgD isotype - IGHV genes")
legend("topleft",legend=c("NP","PNR","PR","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

vusage_sign<-vars

splitop<-strsplit(as.character(vgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))

annotation_col = data.frame(
  sample = factor(subjects),
  time = factor(vgenes_filter$time),
  clin = vgenes_filter$clin)


rownames(annotation_col)<-rownames(vusage_filter2)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   sample= c("sample1" = fill3[1],"sample2" = fill3[2],"sample3"= fill3[3],"sample4"= fill3[4],"sample5"= fill3[5],"sample6"= fill3[6],
                           "sample7"= fill3[7],"sample8"= fill3[8],"sample9"= fill3[9],"sample10"=fill3[10], "sample11"=fill3[11]  ,"sample12"=fill3[12], "sample13"=fill3[13] ,
                             "sample14"=  fill3[14],"sample15"=  fill3[15], "sample16"= fill3[16],"sample17"= fill3[17],"sample18"= fill3[18],"sample19"= fill3[19],"sample20"= fill3[20],
                              "sample21"= fill3[21], "sample22"= fill3[22], "sample23"= fill3[23], "sample24"= fill3[24], "sample25"= fill3[25], "sample26"=fill3[26] , 
                               "sample27"= fill3[27]))


colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_cDNA_IGHD_RF.tiff",res=300,w=3500,h=2500)
pheatmap(t(vusage_sign),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

############
# 2. IGHA
###########
data_cDNA_long_IGHA<-data_cDNA_long[which(data_cDNA_long$isotype=="IGHA"),]
reads_clones_annot_Long_IGHA<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_IGHA>100),]
id<-match(data_cDNA_long_IGHA$specimen_label,reads_clones_annot_Long_IGHA$specimen_id)

data_cDNA_long_IGHA_qc<-data_cDNA_long_IGHA[which(is.na(id)==F),]
data_cDNA_long_IGHA_qc$specimen_label<-factor(data_cDNA_long_IGHA_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_cDNA_long_IGHA_qc$specimen_label,data_cDNA_long_IGHA_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_IGHA$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_IGHA$clin[id.spec],reads_clones_annot_Long_IGHA$time[id.spec],
              reads_clones_annot_Long_IGHA$Sample_id[id.spec],reads_clones_annot_Long_IGHA$subject_id[id.spec])
colnames(vgenes)[60:63]<-c("clin","time","Sample_id","subject_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-4))
for (i in 1:59){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_IGHA$clones_IGHA[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:59]
rownames(vusage)<-vgenes$subject_id


#Filter by time
##First manually convert the time points
vgenes$time<-replace(vgenes$time,vgenes$time==32,24)
vgenes$time<-replace(vgenes$time,vgenes$time==13,24)
vgenes$time<-replace(vgenes$time,vgenes$time==12,6)

##Filter for those that have not time 6 or 24
vgenes_filter<-vgenes[which(vgenes$time==6 | vgenes$time==24),]
vusage_filter<-vusage[which(vgenes$time==6 | vgenes$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.002,0)
###Those genes who are in lesss than 10%
vusage_filter2<-vusage_filter[,which(apply(xx,2,function(x) sum(x==0))<46)]
##47 genes in total

##V-genes found in the gDNA Analysis
# "IGHV4.34" "IGHV4.31" "IGHV3.11" "IGHV3.23" "IGHV3.9"  "IGHV3.33" "IGHV3.21" "IGHV1.46" "IGHV3.7"  "IGHV3.13" "IGHV3.53" "IGHV4.39"
# "IGHV1.24" "IGHV5.a"  "IGHV3.20" "IGHV4.59" "IGHV1.f" 

#########
## Statistical analysis to find significant vgenes
########
set.seed(112233)
vusage_time<-data.frame(cbind(vusage_filter2,vgenes_filter$time))
fit<-VSURF(x = vusage_time, y=factor(vgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vusage_time[,fit$varselect.interp])
set.seed(1000)
rf_output <- randomForest(factor(vgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vgene_cDNA_IGHA_RF.tiff",res=300,w=2000,h=2000)
vgenes_filter$timeplot<-ifelse(vgenes_filter$time==0,1,
                               ifelse(vgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vgenes_filter$clin),cex=vgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "IgA isotype - IGHV genes")
#legend("topleft",legend=c("NP","PNR","PR","t6","t24"), 
#       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
#       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

vusage_sign<-vars

splitop<-strsplit(as.character(vgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))

annotation_col = data.frame(
  sample = factor(subjects),
  time = factor(vgenes_filter$time),
  clin = vgenes_filter$clin)


rownames(annotation_col)<-rownames(vusage_filter2)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   sample= c("sample1" = fill3[1],"sample2" = fill3[2],"sample3"= fill3[3],"sample4"= fill3[4],"sample5"= fill3[5],"sample6"= fill3[6],
                             "sample7"= fill3[7],"sample8"= fill3[8],"sample9"= fill3[9],"sample10"=fill3[10], "sample11"=fill3[11]  ,"sample12"=fill3[12], "sample13"=fill3[13] ,
                             "sample14"=  fill3[14],"sample15"=  fill3[15], "sample16"= fill3[16],"sample17"= fill3[17],"sample18"= fill3[18],"sample19"= fill3[19],"sample20"= fill3[20],
                             "sample21"= fill3[21], "sample22"= fill3[22], "sample23"= fill3[23], "sample24"= fill3[24], "sample25"= fill3[25], "sample26"=fill3[26] , 
                             "sample27"= fill3[27]))


colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_cDNA_IGHA_RF.tiff",res=300,w=3500,h=2500)
pheatmap(t(vusage_sign),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()


############
# 3. IGHM
###########
data_cDNA_long_IGHM<-data_cDNA_long[which(data_cDNA_long$isotype=="IGHM"),]
reads_clones_annot_Long_IGHM<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_IGHM>100),]
id<-match(data_cDNA_long_IGHM$specimen_label,reads_clones_annot_Long_IGHM$specimen_id)

data_cDNA_long_IGHM_qc<-data_cDNA_long_IGHM[which(is.na(id)==F),]
data_cDNA_long_IGHM_qc$specimen_label<-factor(data_cDNA_long_IGHM_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_cDNA_long_IGHM_qc$specimen_label,data_cDNA_long_IGHM_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_IGHM$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_IGHM$clin[id.spec],reads_clones_annot_Long_IGHM$time[id.spec],
              reads_clones_annot_Long_IGHM$Sample_id[id.spec],reads_clones_annot_Long_IGHM$subject_id[id.spec])
colnames(vgenes)[65:68]<-c("clin","time","Sample_id","subject_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-4))
for (i in 1:64){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_IGHM$clones_IGHM[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:64]
rownames(vusage)<-vgenes$subject_id


#Filter by time
##First manually convert the time points
vgenes$time<-replace(vgenes$time,vgenes$time==32,24)
vgenes$time<-replace(vgenes$time,vgenes$time==13,24)
vgenes$time<-replace(vgenes$time,vgenes$time==12,6)

##Filter for those that have not time 6 or 24
vgenes_filter<-vgenes[which(vgenes$time==6 | vgenes$time==24),]
vusage_filter<-vusage[which(vgenes$time==6 | vgenes$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.002,0)
###Those genes who are in lesss than 10%
vusage_filter2<-vusage_filter[,which(apply(xx,2,function(x) sum(x==0))<43)]
##42 genes in total

##V-genes found in the gDNA Analysis
# "IGHV4.34" "IGHV4.31" "IGHV3.11" "IGHV3.23" "IGHV3.9"  "IGHV3.33" "IGHV3.21" "IGHV1.46" "IGHV3.7"  "IGHV3.13" "IGHV3.53" "IGHV4.39"
# "IGHV1.24" "IGHV5.a"  "IGHV3.20" "IGHV4.59" "IGHV1.f" 

#########
## Statistical analysis to find significant vgenes
########
set.seed(112233)
vusage_time<-data.frame(cbind(vusage_filter2,vgenes_filter$time))
fit<-VSURF(x = vusage_time, y=factor(vgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vusage_time[,fit$varselect.interp])
set.seed(1000)
rf_output <- randomForest(factor(vgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vgene_cDNA_IGHM_RF.tiff",res=300,w=2000,h=2000)
vgenes_filter$timeplot<-ifelse(vgenes_filter$time==0,1,
                               ifelse(vgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vgenes_filter$clin),cex=vgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "IgM isotype - IGHV genes")
legend("topleft",legend=c("NP","PNR","PR","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

vusage_sign<-vars

splitop<-strsplit(as.character(vgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))

annotation_col = data.frame(
  sample = factor(subjects),
  time = factor(vgenes_filter$time),
  clin = vgenes_filter$clin)


rownames(annotation_col)<-rownames(vusage_filter2)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   sample= c("sample1" = fill3[1],"sample2" = fill3[2],"sample3"= fill3[3],"sample4"= fill3[4],"sample5"= fill3[5],"sample6"= fill3[6],
                             "sample7"= fill3[7],"sample8"= fill3[8],"sample9"= fill3[9],"sample10"=fill3[10], "sample11"=fill3[11]  ,"sample12"=fill3[12], "sample13"=fill3[13] ,
                             "sample14"=  fill3[14],"sample15"=  fill3[15], "sample16"= fill3[16],"sample17"= fill3[17],"sample18"= fill3[18],"sample19"= fill3[19],"sample20"= fill3[20],
                             "sample21"= fill3[21], "sample22"= fill3[22], "sample23"= fill3[23], "sample24"= fill3[24], "sample25"= fill3[25], "sample26"=fill3[26] , 
                             "sample27"= fill3[27]))


colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_cDNA_IGHM_RF.tiff",res=300,w=3500,h=2500)
pheatmap(t(vusage_sign),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

############
# 4. IGHG
###########
data_cDNA_long_IGHG<-data_cDNA_long[which(data_cDNA_long$isotype=="IGHG"),]
reads_clones_annot_Long_IGHG<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_IGHG>100),]
id<-match(data_cDNA_long_IGHG$specimen_label,reads_clones_annot_Long_IGHG$specimen_id)

data_cDNA_long_IGHG_qc<-data_cDNA_long_IGHG[which(is.na(id)==F),]
data_cDNA_long_IGHG_qc$specimen_label<-factor(data_cDNA_long_IGHG_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_cDNA_long_IGHG_qc$specimen_label,data_cDNA_long_IGHG_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_IGHG$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_IGHG$clin[id.spec],reads_clones_annot_Long_IGHG$time[id.spec],
              reads_clones_annot_Long_IGHG$Sample_id[id.spec],reads_clones_annot_Long_IGHG$subject_id[id.spec])
colnames(vgenes)[60:63]<-c("clin","time","Sample_id","subject_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-4))
for (i in 1:59){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_IGHG$clones_IGHG[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:59]
rownames(vusage)<-vgenes$subject_id


#Filter by time
##First manually convert the time points
vgenes$time<-replace(vgenes$time,vgenes$time==32,24)
vgenes$time<-replace(vgenes$time,vgenes$time==13,24)
vgenes$time<-replace(vgenes$time,vgenes$time==12,6)

##Filter for those that have not time 6 or 24
vgenes_filter<-vgenes[which(vgenes$time==6 | vgenes$time==24),]
vusage_filter<-vusage[which(vgenes$time==6 | vgenes$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.002,0)
###Those genes who are in lesss than 10%
vusage_filter2<-vusage_filter[,which(apply(xx,2,function(x) sum(x==0))<46)]
##47 genes in total

##V-genes found in the gDNA Analysis
# "IGHV4.34" "IGHV4.31" "IGHV3.11" "IGHV3.23" "IGHV3.9"  "IGHV3.33" "IGHV3.21" "IGHV1.46" "IGHV3.7"  "IGHV3.13" "IGHV3.53" "IGHV4.39"
# "IGHV1.24" "IGHV5.a"  "IGHV3.20" "IGHV4.59" "IGHV1.f" 

#########
## Statistical analysis to find significant vgenes
########
set.seed(112233)
vusage_time<-data.frame(cbind(vusage_filter2,vgenes_filter$time))
fit<-VSURF(x = vusage_time, y=factor(vgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vusage_time[,fit$varselect.interp])
set.seed(1000)
rf_output <- randomForest(factor(vgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vgene_cDNA_IGHG_RF.tiff",res=300,w=2000,h=2000)
vgenes_filter$timeplot<-ifelse(vgenes_filter$time==0,1,
                               ifelse(vgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vgenes_filter$clin),cex=vgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "IgG isotype - IGHV genes")
legend("bottomleft",legend=c("NP","PNR","PR","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

vusage_sign<-vars

splitop<-strsplit(as.character(vgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))

annotation_col = data.frame(
  sample = factor(subjects),
  time = factor(vgenes_filter$time),
  clin = vgenes_filter$clin)


rownames(annotation_col)<-rownames(vusage_filter2)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   sample= c("sample1" = fill3[1],"sample2" = fill3[2],"sample3"= fill3[3],"sample4"= fill3[4],"sample5"= fill3[5],"sample6"= fill3[6],
                             "sample7"= fill3[7],"sample8"= fill3[8],"sample9"= fill3[9],"sample10"=fill3[10], "sample11"=fill3[11]  ,"sample12"=fill3[12], "sample13"=fill3[13] ,
                             "sample14"=  fill3[14],"sample15"=  fill3[15], "sample16"= fill3[16],"sample17"= fill3[17],"sample18"= fill3[18],"sample19"= fill3[19],"sample20"= fill3[20],
                             "sample21"= fill3[21], "sample22"= fill3[22], "sample23"= fill3[23], "sample24"= fill3[24], "sample25"= fill3[25], "sample26"=fill3[26] , 
                             "sample27"= fill3[27]))


colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_cDNA_IGHG_RF.tiff",res=300,w=3500,h=2500)
pheatmap(t(vusage_sign),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()