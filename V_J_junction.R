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

setwd("/Users/Pinedasans/VDJ/SummaryResults/Vgene/")
load("/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")

############
### gDNA ###
############
data_gDNA<-data_merge[which(data_merge$amplification_template=="gDNA"),]
### longitudinal
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]

#For this analysis we are putting a cut-off on clones because v-genes can be biased at low clonality 
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_gDNA>100),]

id<-match(data_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_gDNA_long_qc$specimen_label,data_gDNA_long_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_qc$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time[id.spec],
              reads_clones_annot_Long_qc$Sample_id[id.spec],reads_clones_annot_Long_qc$subject_id[id.spec])
colnames(vgenes)[65:68]<-c("clin","time","Sample_id","subject_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-4))
for (i in 1:64){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:64]
rownames(vusage)<-vgenes$subject_id

##Delete sample 8
vgenes2<-vgenes[which(vgenes$subject_id!="sample8_6" & vgenes$subject_id!="sample8_24"),]
vusage2<-vusage[which(vgenes$subject_id!="sample8_6" & vgenes$subject_id!="sample8_24"),]

#Replace some time points
vgenes2$time<-replace(vgenes2$time,vgenes2$time==13,24)
vgenes2$time<-replace(vgenes2$time,vgenes2$time==12,6)
vgenes2$time<-replace(vgenes2$time,vgenes2$time==32,24)

#to delete a sample with time==9 and ==2
vgenes_filter<-vgenes2[which(vgenes2$time==0 | vgenes2$time==6 | vgenes2$time==24),]
vusage_filter<-vusage2[which(vgenes2$time==0 | vgenes2$time==6 | vgenes2$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.002,0)
###Those who are in lesss than 10%
vusage_filter<-vusage_filter[,which(apply(xx,2,function(x) sum(x==0))<=57)]
##47 genes in total

#########
## Statistical analysis to find significant vgenes
########

##ANOVA by time points
p.PNR<-matrix(NA,4,64)
p.PR<-matrix(NA,4,64)
for(i in 1:64){
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
###No signifcances by time point

###Longitudinal analysis
p<-NULL
for(i in 1:64){
  fm_null <- lmer(vusage_filter[,i] ~ vgenes_filter$clin + vgenes_filter$time + (vgenes_filter$time | vgenes_filter$Sample_id),REML = F)
  fm_full <- lmer(vusage_filter[,i] ~  vgenes_filter$clin*vgenes_filter$time + (vgenes_filter$time | vgenes_filter$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##
} 
##Non significant results across time points

###Rnadom Forest
set.seed(112233)
vusage_time<-data.frame(cbind(vusage_filter,vgenes_filter$time))
fit<-VSURF(x = vusage_time, y=factor(vgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vusage_time[,fit$varselect.interp])

###13 and 43 has very low expression, deleted
#vars<-vars[,which(colnames(vars)!="IGHV3.13" & colnames(vars)!="IGHV3.43")]

set.seed(1000)
rf_output <- randomForest(factor(vgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vgene_RF_gDNA.tiff",res=300,w=2000,h=2000)
vgenes_filter$timeplot<-ifelse(vgenes_filter$time==0,1,
                               ifelse(vgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vgenes_filter$clin),cex=vgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "IGHV genes")
legend("topleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()


##Heatmap with significant results
vusage_sign<-vars

splitop<-strsplit(as.character(vgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))
individuals<-c("individual1" , "individual1" , "individual1",  "individual2" , "individual2"  ,"individual3" , "individual3" , "individual4"  ,"individual4" , "individual4" , "individual5" 
               ,"individual6" , "individual6" ,"individual7" , "individual8" , "individual9",  "individual9" , "individual9" , "individual10" ,"individual10", "individual10", "individual11" 
               ,"individual11", "individual12" ,"individual12" ,"individual12", "individual13", "individual13" ,"individual13" ,"individual14", "individual14" ,"individual16" ,"individual16" 
               ,"individual16" ,"individual17" ,"individual17", "individual18", "individual18", "individual19","individual19" ,"individual19" ,"individual20", "individual20", "individual20" 
               ,"individual21" ,"individual21", "individual21", "individual22", "individual22", "individual22", "individual23" ,"individual23", "individual23", "individual24", "individual24" 
               ,"individual24", "individual25" ,"individual25" ,"individual25" ,"individual26", "individual26" ,"individual26" ,"individual27")

annotation_col = data.frame(
  individual = factor(individuals),
  time = as.factor(vgenes_filter$time),
  clin = vgenes_filter$clin)


rownames(annotation_col)<-rownames(vusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- 26
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   individual= c("individual1" = fill3[1],"individual2" = fill3[2],"individual3"= fill3[3],"individual4"= fill3[4],"individual5"= fill3[5],"individual6"= fill3[6],
                           "individual7"= fill3[7],"individual8"= fill3[8],"individual9"= fill3[9],"individual10"=fill3[10], "individual11"=fill3[11]  ,"individual12"=fill3[12], "individual13"=fill3[13] ,
                             "individual14"=  fill3[14], "individual16"= fill3[15],"individual17"= fill3[16],"individual18"= fill3[17],"individual19"= fill3[18],"individual20"= fill3[19],
                              "individual21"= fill3[20], "individual22"= fill3[21], "individual23"= fill3[22], "individual24"= fill3[23], "individual25"= fill3[24], "individual26"=fill3[25] , 
                               "individual27"= fill3[26]))


colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_byclones_RF_gDNA.tiff",res=300,w=3500,h=2500)
pheatmap(t(vusage_sign),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

#################################
####### V  J junction #########
###############################
data_qc_gDNA_long_qc$vjgenes<- paste(data_qc_gDNA_long_qc$v_gene,data_qc_gDNA_long_qc$j_gene,sep="_")
vjgenes<-as.data.frame(unclass(table(data_qc_gDNA_long_qc$specimen_label,data_qc_gDNA_long_qc$vjgenes)))
id.spec<-match(rownames(vjgenes),reads_clones_annot_Long_qc$specimen_id)
vjgenes<-cbind(vjgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec],reads_clones_annot_Long_qc$subject_id.1[id.spec])
colnames(vjgenes)[326:329]<-c("clin","time","Sample_id","subject_id")

vjusage<-matrix(NA,nrow(vjgenes),(ncol(vjgenes)-4))
for (i in 1:325){
  vjusage[,i]<-vjgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(vjusage)<-colnames(vjgenes)[1:325]
rownames(vjusage)<-vjgenes$subject_id

##Delete sample 8
vjgenes2<-vjgenes[which(vjgenes$subject_id!="sample8_6" & vjgenes$subject_id!="sample8_24"),]
vjusage2<-vjusage[which(vjgenes$subject_id!="sample8_6" & vjgenes$subject_id!="sample8_24"),]

#to delete a sample with time==9 and ==2
vjgenes_filter<-vjgenes2[which(vjgenes2$time==0 | vjgenes2$time==6 | vjgenes2$time==24),]
vjusage_filter<-vjusage2[which(vjgenes2$time==0 | vjgenes2$time==6 | vjgenes2$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vjusage_filter,vjusage_filter<0.002,0)
#xx[,which(apply(xx,2,function(x) sum(x==0))<=62)]
###Those who are in lesss than 10%
vjusage_filter<-vjusage_filter[,which(apply(xx,2,function(x) sum(x==0))<=58)]
##213 vj junction 



#########
## Statistical analysis to find significant jgenes
########
p<-NULL
for(i in 1:284){
  fm_null <- lmer(vjusage_filter[,i] ~ vjgenes2$clin + vjgenes2$time + (vjgenes2$time | vjgenes2$Sample_id),REML = F)
  fm_full <- lmer(vjusage_filter[,i] ~  vjgenes2$clin*vjgenes2$time + (vjgenes2$time | vjgenes2$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##non
  
} 

set.seed(1234)
vjusage_time<-cbind(vjusage_filter,vjgenes_filter$time)
fit<-VSURF(x = vjusage_time, y=factor(vjgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vjusage_time[,fit$varselect.interp])
set.seed(400)
rf_output <- randomForest(factor(vjgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vj_gene_RF.tiff",res=300,w=2500,h=2500)
vjgenes_filter$timeplot<-ifelse(vjgenes_filter$time==0,1,
                                 ifelse(vjgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vjgenes_filter$clin),cex=vjgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "V-J junction")
legend("topleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()


##Heatmap
vjusage_sign<-vars

splitop<-strsplit(as.character(vjgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))
individuals<-c("individual1" , "individual1" , "individual1",  "individual2" , "individual2"  ,"individual3" , "individual3" , "individual4"  ,"individual4" , "individual4" , "individual5" 
               ,"individual6" , "individual6" ,"individual7" , "individual8" , "individual9",  "individual9" , "individual9" , "individual10" ,"individual10", "individual10", "individual11" 
               ,"individual11", "individual12" ,"individual12" ,"individual12", "individual13", "individual13" ,"individual13" ,"individual14", "individual14" ,"individual16" ,"individual16" 
               ,"individual16" ,"individual17" ,"individual17", "individual18", "individual18", "individual19","individual19" ,"individual19" ,"individual20", "individual20", "individual20" 
               ,"individual21" ,"individual21", "individual21", "individual22", "individual22", "individual22", "individual23" ,"individual23", "individual23", "individual24", "individual24" 
               ,"individual24", "individual25" ,"individual25" ,"individual25" ,"individual26", "individual26" ,"individual26" ,"individual27")

annotation_col = data.frame(
  individual = factor(individuals),
  time = as.factor(vjgenes_filter$time),
  clin = vjgenes_filter$clin)


rownames(annotation_col)<-rownames(vjusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- 26
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   individual= c("individual1" = fill3[1],"individual2" = fill3[2],"individual3"= fill3[3],"individual4"= fill3[4],"individual5"= fill3[5],"individual6"= fill3[6],
                                 "individual7"= fill3[7],"individual8"= fill3[8],"individual9"= fill3[9],"individual10"=fill3[10], "individual11"=fill3[11]  ,"individual12"=fill3[12], "individual13"=fill3[13] ,
                                 "individual14"=  fill3[14], "individual16"= fill3[15],"individual17"= fill3[16],"individual18"= fill3[17],"individual19"= fill3[18],"individual20"= fill3[19],
                                 "individual21"= fill3[20], "individual22"= fill3[21], "individual23"= fill3[22], "individual24"= fill3[23], "individual25"= fill3[24], "individual26"=fill3[25] , 
                                 "individual27"= fill3[26]))

colfunc<-colorRampPalette(c("white","purple4"))                
tiff("heatmap_vdjgene_RF.tiff",res=300,w=3500,h=2500)
pheatmap(t(vjusage_sign),cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()



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






