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
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")

############
### cDNA ###
############
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
### longitudinal
data_qc_cDNA_long<-data_qc_cDNA[which(data_qc_cDNA$clin!="AR" & data_qc_cDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]

#For this analysis we are putting a cut-off on clones because v-genes can be biased at low clonality 
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_cDNA>1000),]

id<-match(data_qc_cDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_cDNA_long_qc<-data_qc_cDNA_long[which(is.na(id)==F),]
data_qc_cDNA_long_qc$specimen_label<-factor(data_qc_cDNA_long_qc$specimen_label)

vgenes<-as.data.frame(unclass(table(data_qc_cDNA_long_qc$specimen_label,data_qc_cDNA_long_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_qc$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec],reads_clones_annot_Long_qc$subject_id.1[id.spec])
colnames(vgenes)[70:73]<-c("clin","time","Sample_id","subject_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-4))
for (i in 1:69){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_qc$clones_cDNA[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:69]
rownames(vusage)<-vgenes$subject_id


#to delete samples that are not time0,6,24
vgenes_filter<-vgenes[which(vgenes$time==0 | vgenes$time==6 | vgenes$time==24),]
vusage_filter<-vusage[which(vgenes$time==0 | vgenes$time==6 | vgenes$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage_filter,vusage_filter<0.002,0)
###Those who are in lesss than 10%
vusage_filter<-vusage_filter[,which(apply(vusage_filter,2,function(x) sum(x==0))<=46)]
##61 genes in total


#########
## Statistical analysis to find significant vgenes
########
set.seed(112233)
vusage_time<-data.frame(cbind(vusage_filter,vgenes_filter$time))
fit<-VSURF(x = vusage_time, y=factor(vgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vusage_time[,fit$varselect.interp])
set.seed(1000)
rf_output <- randomForest(factor(vgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vgene_RF.tiff",res=300,w=2000,h=2000)
vgenes_filter$timeplot<-ifelse(vgenes_filter$time==0,1,
                               ifelse(vgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vgenes_filter$clin),cex=vgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "IGHV genes")
legend("topright",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

vusage_sign<-vars

splitop<-strsplit(as.character(vgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))

annotation_col = data.frame(
  sample = factor(subjects),
  time = as.factor(vgenes_filter$time),
  clin = vgenes_filter$clin)


rownames(annotation_col)<-rownames(vusage_filter)
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
tiff("heatmap_vgene_byclones_RF.tiff",res=300,w=3500,h=2500)
pheatmap(t(vusage_sign),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

#################################
####### V  J junction #########
###############################
data_qc_cDNA_long_qc$vjgenes<- paste(data_qc_cDNA_long_qc$v_gene,data_qc_cDNA_long_qc$j_gene,sep="_")
vjgenes<-as.data.frame(unclass(table(data_qc_cDNA_long_qc$specimen_label,data_qc_cDNA_long_qc$vjgenes)))
id.spec<-match(rownames(vjgenes),reads_clones_annot_Long_qc$specimen_id)
vjgenes<-cbind(vjgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec],reads_clones_annot_Long_qc$subject_id.1[id.spec])
colnames(vjgenes)[353:356]<-c("clin","time","Sample_id","subject_id")

vjusage<-matrix(NA,nrow(vjgenes),(ncol(vjgenes)-4))
for (i in 1:352){
  vjusage[,i]<-vjgenes[,i]/reads_clones_annot_Long_qc$clones_cDNA[id.spec]
}
colnames(vjusage)<-colnames(vjgenes)[1:352]
rownames(vjusage)<-vjgenes$subject_id


#to delete a sample with time==9 and ==2
vjgenes_filter<-vjgenes[which(vjgenes$time==0 | vjgenes$time==6 | vjgenes$time==24),]
vjusage_filter<-vjusage[which(vjgenes$time==0 | vjgenes$time==6 | vjgenes$time==24),]

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vjusage_filter,vjusage_filter<0.002,0)
###Those who are in lesss than 10%
vjusage_filter<-vjusage_filter[,which(apply(xx,2,function(x) sum(x==0))<=46)]
##205 vj junction 



#########
## Statistical analysis to find significant jgenes
########

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

vjusage_sign<-vars

splitop<-strsplit(as.character(vjgenes_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))

annotation_col = data.frame(
  sample = factor(subjects),
  time = as.factor(vjgenes_filter$time),
  clin = vjgenes_filter$clin)


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
names(grid.col)<-c(as.character(unique(v)),as.character(unique(j)))


tiff(paste("v_j_junction_cDNA_long_PR.tiff"),res=300,h=1800,w=1800)

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

names(grid.col)<-c(as.character(unique(v)),as.character(unique(j)))

tiff(paste("v_j_junction_cDNA_long_PNR.tiff"),res=300,h=1800,w=1800)

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
names(grid.col)<-c(as.character(unique(v)),as.character(unique(j)))

tiff(paste("v_j_junction_cDNA_long_NP.tiff"),res=300,h=1800,w=1800)

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






