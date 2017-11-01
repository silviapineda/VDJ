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

setwd("/Users/Pinedasans/VDJ/Results/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")

############
### gDNA ###
############
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
### longitudinal
data_qc_gDNA_long<-data_qc_gDNA[which(data_qc_gDNA$clin!="AR" & data_qc_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_gDNA>100),]

id<-match(data_qc_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_gDNA_long_qc<-data_qc_gDNA_long[which(is.na(id)==F),]
data_qc_gDNA_long_qc$specimen_label<-factor(data_qc_gDNA_long_qc$specimen_label)


####### V genes #########
vgenes<-as.data.frame(unclass(table(data_qc_gDNA_long_qc$specimen_label,data_qc_gDNA_long_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_Long_qc$specimen_id)
vgenes<-cbind(vgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec])
colnames(vgenes)[65:67]<-c("clin","time","Sample_id")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-3))
for (i in 1:64){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:64]
rownames(vusage)<-rownames(vgenes)

vgenes_filter<-vgenes[which(vgenes$time==0 | vgenes$time==6 | vgenes$time==24),]
vusage_filter<-vusage[which(vgenes$time==0 | vgenes$time==6 | vgenes$time==24),]

###At least 5% of the samples has expression
vusage_filter2<-vusage_filter[,which(apply(vusage_filter,2,function(x) sum(x==0))<=64)]

#normalized = (x-min(x))/(max(x)-min(x))
#vusage_norm<-apply(vusage_filter,1,function(x) (x-min(x))/(max(x)-min(x)))


#### heatmap
#x <- t(scale(vusage)) #scale genes by columns, samples by rows

annotation_col = data.frame(
  sample = factor(vgenes_filter$Sample_id),
  time = as.factor(vgenes_filter$time),
  clin = vgenes_filter$clin)
  
rownames(annotation_col)<-rownames(vusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
fill3=rainbow(26)
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),
                   sample= c("1005" = fill3[1],"1006" = fill3[2],"12001"= fill3[3],"12004"= fill3[4],"12008"= fill3[5],"38010"= fill3[6],
                          "38011"= fill3[7],"38012"= fill3[8],"38020"= fill3[9],"38021"=fill3[10], "41006"=fill3[11]  ,"89002"=fill3[12], "300003"=fill3[13] ,
                            "300005"=  fill3[14], "300006"= fill3[15],"300013"= fill3[16],"300015"= fill3[17],"301002"= fill3[18],"301003"= fill3[19],
                          "303001"= fill3[20], "303003"= fill3[21], "303005"= fill3[22], "303006"= fill3[23], "303007"= fill3[24], "303010"=fill3[25] , 
                           "303013"= fill3[26]))
tiff("heatmap_vgene_byclones.tiff",res=300,w=4000,h=2500)
pheatmap(t(vusage_filter2),cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Reds"),border_color=F)
dev.off()


#########
## Statistical analysis to find significant vgenes
########

##ANOVA by time points
p.PNR<-matrix(NA,4,62)
p.PR<-matrix(NA,4,62)
for(i in 1:62){
  fit0<-summary(glm(vusage_filter[which(vgenes_filter$time==0),i] ~ vgenes_filter$clin[which(vgenes_filter$time==0)]))
  p.PNR[1,i]<-fit0$coefficients[2,4]
  p.PR[1,i]<-fit0$coefficients[3,4]
  fit6<-summary(glm(vusage_filter[which(vgenes_filter$time==6),i] ~ vgenes_filter$clin[which(vgenes_filter$time==6)]))
  p.PNR[2,i]<-fit6$coefficients[2,4]
  p.PR[2,i]<-fit6$coefficients[3,4]
  fit24<-summary(glm(vusage_filter[which(vgenes_filter$time==24),i] ~ vgenes_filter$clin[which(vgenes_filter$time==24)]))
  p.PNR[3,i]<-fit24$coefficients[2,4]
  p.PR[3,i]<-fit24$coefficients[3,4]
  fit<-summary(glm(vusage_filter[,i] ~ vgenes_filter$clin))
  p.PNR[4,i]<-fit$coefficients[2,4]
  p.PR[4,i]<-fit$coefficients[3,4]
}
##PNR not statistically significance
colnames(p.PNR)<-colnames(vusage_filter2)
colnames(p.PR)<-colnames(vusage_filter2)
p.PR[,c("IGHV1-18","IGHV1-24","IGHV1-58","IGHV1-8","IGHV2-26","IGHV2-70","IGHV3-11","IGHV3-13","IGHV3-15","IGHV3-23","IGHV3-30","IGHV3-33",
        "IGHV3-43","IGHV3-47","IGHV3-53","IGHV3-7","IGHV3-9","IGHV4-34","IGHV4-39")]

vusage_sign<-data.frame(vusage_filter2[,match(c("IGHV1-18","IGHV1-24","IGHV1-58","IGHV1-8","IGHV2-26","IGHV2-70","IGHV3-11","IGHV3-13","IGHV3-15","IGHV3-23","IGHV3-30","IGHV3-33",
        "IGHV3-43","IGHV3-47","IGHV3-53","IGHV3-7","IGHV3-9","IGHV4-34","IGHV4-39"),colnames(vusage))])

vusage_sign_0<-vusage_sign[which(vgenes_filter$time==0),]
vgenes_filter_0<-vgenes_filter[which(vgenes_filter$time==0),]
vusage_NP_0<-apply(vusage_sign_0[which(vgenes_filter_0$clin=="NP"),],2,function (x) mean(x))
vusage_PNR_0<-apply(vusage_sign_0[which(vgenes_filter_0$clin=="PNR"),],2,function (x) mean(x))
vusage_PR_0<-apply(vusage_sign_0[which(vgenes_filter_0$clin=="PR"),],2,function (x) mean(x))

vusage_sign_6<-vusage_sign[which(vgenes_filter$time==6),]
vgenes_filter_6<-vgenes_filter[which(vgenes_filter$time==6),]
vusage_NP_6<-apply(vusage_sign_6[which(vgenes_filter_6$clin=="NP"),],2,function (x) mean(x))
vusage_PNR_6<-apply(vusage_sign_6[which(vgenes_filter_6$clin=="PNR"),],2,function (x) mean(x))
vusage_PR_6<-apply(vusage_sign_6[which(vgenes_filter_6$clin=="PR"),],2,function (x) mean(x))

vusage_sign_24<-vusage_sign[which(vgenes_filter$time==24),]
vgenes_filter_24<-vgenes_filter[which(vgenes_filter$time==24),]
vusage_NP_24<-apply(vusage_sign_24[which(vgenes_filter_24$clin=="NP"),],2,function (x) mean(x))
vusage_PNR_24<-apply(vusage_sign_24[which(vgenes_filter_24$clin=="PNR"),],2,function (x) mean(x))
vusage_PR_24<-apply(vusage_sign_24[which(vgenes_filter_24$clin=="PR"),],2,function (x) mean(x))

vusage_NP<-apply(vusage_sign[which(vgenes_filter$clin=="NP"),],2,function (x) mean(x))
vusage_PNR<-apply(vusage_sign[which(vgenes_filter$clin=="PNR"),],2,function (x) mean(x))
vusage_PR<-apply(vusage_sign[which(vgenes_filter$clin=="PR"),],2,function (x) mean(x))

# Simple Pie Chart
tiff("PieChart_VgenesSign2.tiff",res=300,w=1850,h=1850)
fill<-c("#e6194b",	"#3cb44b",	"#ffe119",	"#0082c8",	"#f58231",	"#911eb4",	"#46f0f0",	"#f032e6",	"#d2f53c",
        "#fabebe",	"#008080",	"#e6beff",	"#aa6e28",	"#fffac8",	"#800000",	"#aaffc3",	"#808000",	"#ffd8b1",
        "#000080",	"#808080",	"#FFFFFF",	"#000000")
par(mfrow=c(3,3))
pie(vusage_NP_0, labels = NA,radius=1,col=fill,clockwise = T)
pie(vusage_PNR_0, labels = NA, radius=1,col=fill,clockwise = T)
pie(vusage_PR_0, labels = NA, radius=1,col=fill,clockwise = T)

pie(vusage_NP_6, labels = NA,radius=1,col=fill,clockwise = T)
pie(vusage_PNR_6, labels = NA,radius=1,col=fill,clockwise = T)
pie(vusage_PR_6, labels = NA,radius=1,col=fill,clockwise = T)

pie(vusage_NP_24, labels = NA, radius=1,col=fill,clockwise = T)
pie(vusage_PNR_24, labels = NA, radius=1,col=fill,clockwise = T)
pie(vusage_PR_24, labels = NA, radius=1,col=fill,clockwise = T)

legend(x=-13.7, y=8.5, legend=names(vusage_NP), fill=fill,xpd=NA, bty="n")
dev.off()

tiff("PieChart_VgenesSignTotal.tiff",res=300,w=2000,h=2000)
par(mfrow=c(1,3))
pie(vusage_NP, labels = NA, radius=1,col=fill,clockwise = T)
pie(vusage_PNR, labels = NA, radius=1,col=fill,clockwise = T)
pie(vusage_PR, labels = NA, radius=1,col=fill,clockwise = T)

#legend("bottom", legend=names(vusage_NP), fill=fill,xpd=NA, bty="n")
dev.off()

p<-NULL
for(i in 1:62){
  fm_null <- lmer(vusage_filter2[,i] ~ vgenes_filter$clin + vgenes_filter$time + (vgenes_filter$time | vgenes_filter$Sample_id),REML = F)
  fm_full <- lmer(vusage_filter2[,i] ~  vgenes_filter$clin*vgenes_filter$time + (vgenes_filter$time | vgenes_filter$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##
} 

set.seed(4000)
vusage_time<-data.frame(cbind(vusage_filter,vgenes_filter$time))
fit<-VSURF(x = vusage_time, y=factor(vgenes_filter$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vusage_time[,fit$varselect.interp])
set.seed(1000)
rf_output <- randomForest(factor(vgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vgene_RF.tiff",res=300,w=2000,h=2000)
vgenes_filter$timeplot<-ifelse(vgenes_filter$time==0,1,
                               ifelse(vgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vgenes_filter$clin),cex=vgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "RF selection")
legend("topleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

vusage_sign<-vars
#vusage_norm<-apply(vusage_sign,1,function(x) (x-min(x))/(max(x)-min(x)))

annotation_col = data.frame(
  sample = factor(vgenes_filter$Sample_id),
  time = as.factor(vgenes_filter$time),
  clin = vgenes_filter$clin)

rownames(annotation_col)<-rownames(vusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
fill3=rainbow(26)
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),
                   sample= c("1005" = fill3[1],"1006" = fill3[2],"12001"= fill3[3],"12004"= fill3[4],"12008"= fill3[5],"38010"= fill3[6],
                             "38011"= fill3[7],"38012"= fill3[8],"38020"= fill3[9],"38021"=fill3[10], "41006"=fill3[11]  ,"89002"=fill3[12], "300003"=fill3[13] ,
                             "300005"=  fill3[14], "300006"= fill3[15],"300013"= fill3[16],"300015"= fill3[17],"301002"= fill3[18],"301003"= fill3[19],
                             "303001"= fill3[20], "303003"= fill3[21], "303005"= fill3[22], "303006"= fill3[23], "303007"= fill3[24], "303010"=fill3[25] , 
                             "303013"= fill3[26]))
tiff("heatmap_vgene_byclones_RF.tiff",res=300,w=3000,h=2000)
pheatmap(t(vusage_sign),cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Reds"),border_color=F)
dev.off()


####### V D J junction #########
data_qc_gDNA_long_qc$vdjgenes<- paste(data_qc_gDNA_long_qc$v_gene,data_qc_gDNA_long_qc$d_gene, data_qc_gDNA_long_qc$j_gene,sep="_")
vdjgenes<-as.data.frame(unclass(table(data_qc_gDNA_long_qc$specimen_label,data_qc_gDNA_long_qc$vdjgenes)))
id.spec<-match(rownames(vdjgenes),reads_clones_annot_Long_qc$specimen_id)
vdjgenes<-cbind(vdjgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec])
colnames(vdjgenes)[6097:6099]<-c("clin","time","Sample_id")

vdjusage<-matrix(NA,nrow(vdjgenes),(ncol(vdjgenes)-3))
for (i in 1:6096){
  vdjusage[,i]<-vdjgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(vdjusage)<-colnames(vdjgenes)[1:6096]
rownames(vdjusage)<-rownames(vdjgenes)

vdjgenes_filter<-vdjgenes[which(vdjgenes$time==0 | vdjgenes$time==6 | vdjgenes$time==24),]
vdjusage_filter<-vdjusage[which(vdjgenes$time==0 | vdjgenes$time==6 | vdjgenes$time==24),]

###At least 5% of the samples has expression
vdjusage_filter2<-vdjusage_filter[,which(apply(vdjusage_filter,2,function(x) sum(x==0))<=64)]

#### heatmap
annotation_col = data.frame(
  sample = factor(vdjgenes_filter$Sample_id),
  time = as.factor(vdjgenes_filter$time),
  clin = vdjgenes_filter$clin)

rownames(annotation_col)<-rownames(vdjusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
fill3=rainbow(26)
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),
                   sample= c("1005" = fill3[1],"1006" = fill3[2],"12001"= fill3[3],"12004"= fill3[4],"12008"= fill3[5],"38010"= fill3[6],
                             "38011"= fill3[7],"38012"= fill3[8],"38020"= fill3[9],"38021"=fill3[10], "41006"=fill3[11]  ,"89002"=fill3[12], "300003"=fill3[13] ,
                             "300005"=  fill3[14], "300006"= fill3[15],"300013"= fill3[16],"300015"= fill3[17],"301002"= fill3[18],"301003"= fill3[19],
                             "303001"= fill3[20], "303003"= fill3[21], "303005"= fill3[22], "303006"= fill3[23], "303007"= fill3[24], "303010"=fill3[25] , 
                             "303013"= fill3[26]))
tiff("heatmap_vdjgene.tiff",res=300,w=4000,h=2500)
pheatmap(t(vdjusage_filter2),cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Reds"),border_color=F)
dev.off()


#########
## Statistical analysis to find significant jgenes
########
p<-NULL
for(i in 1:6096){
  fm_null <- lmer(v_d_j_usage[,i] ~ v_d_j_genes$clin + v_d_j_genes$time + (v_d_j_genes$time | v_d_j_genes$Sample_id),REML = F)
  fm_full <- lmer(v_d_j_usage[,i] ~  v_d_j_genes$clin*v_d_j_genes$time + (v_d_j_genes$time | v_d_j_genes$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##non
  
} 

set.seed(4000)
vdjusage_time<-cbind(vdjusage_filter2,vdjgenes_filter$time)
save(vdjusage_time,vdjgenes_filter,file="vdj_reads.Rdata")
fit<-VSURF(x = vdjusage_time, y=factor(vdjgenes_filter$clin),parallel = TRUE,ncores=4)
#save(fit,file="RF_v_d_j_usage.Rdata")
load("RF_vdjusage.Rdata")
vars<-data.frame(vdjusage_time[,fit$varselect.interp])
set.seed(400)
rf_output <- randomForest(factor(vdjgenes_filter$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_vdj_gene_RF.tiff",res=300,w=2000,h=2000)
vdjgenes_filter$timeplot<-ifelse(vdjgenes_filter$time==0,1,
                               ifelse(vdjgenes_filter$time==6,1.5,2))
MDSplot(rf_output, factor(vdjgenes_filter$clin),cex=vdjgenes_filter$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "RF selection")
legend("bottomleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

vdjusage_sign<-vars
vdjusage_norm<-apply(vdjusage_sign,1,function(x) (x-min(x))/(max(x)-min(x)))

annotation_col = data.frame(
  sample = factor(vdjgenes_filter$Sample_id),
  time = as.factor(vdjgenes_filter$time),
  clin = vdjgenes_filter$clin)

rownames(annotation_col)<-rownames(vdjusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
fill3=rainbow(26)
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),
                   sample= c("1005" = fill3[1],"1006" = fill3[2],"12001"= fill3[3],"12004"= fill3[4],"12008"= fill3[5],"38010"= fill3[6],
                             "38011"= fill3[7],"38012"= fill3[8],"38020"= fill3[9],"38021"=fill3[10], "41006"=fill3[11]  ,"89002"=fill3[12], "300003"=fill3[13] ,
                             "300005"=  fill3[14], "300006"= fill3[15],"300013"= fill3[16],"300015"= fill3[17],"301002"= fill3[18],"301003"= fill3[19],
                             "303001"= fill3[20], "303003"= fill3[21], "303005"= fill3[22], "303006"= fill3[23], "303007"= fill3[24], "303010"=fill3[25] , 
                             "303013"= fill3[26]))
tiff("heatmap_vdjgene_byreads_RF_norm.tiff",res=300,w=3000,h=2000)
pheatmap(vdjusage_norm,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Oranges"),border_color=F)
dev.off()


#####Run the circos plot for the significant vdj junctions
id_specimen<-rownames(vdjusage_sign)
for(i in 1:length(id_specimen)){
  print(i)
  specimen<-t(vdjusage_sign[i,])
  specimen2<-specimen[which(specimen[,1]!=0),]
  if(length(specimen2)==1){
    names(specimen2)=rownames(specimen)[which(specimen!=0)]
  }
  splitop<-strsplit(names(specimen2),"_")
  v<-unlist(lapply(splitop, `[[`, 1))
  d<-unlist(lapply(splitop, `[[`, 2))
  j<-unlist(lapply(splitop, `[[`, 3))
  vd_usage<-data.frame(cbind(v,d,as.numeric(specimen2)))
  vj_usage<-data.frame(cbind(v,j,as.numeric(specimen2)))
  dj_usage<-data.frame(cbind(d,j,as.numeric(specimen2)))
  colnames(vd_usage)<-c("link1","link2","count")
  colnames(vj_usage)<-c("link1","link2","count")
  colnames(dj_usage)<-c("link1","link2","count")
  vdj_usage_df<-rbind(vd_usage,vj_usage,dj_usage)
  vdj_usage_df$count<-as.numeric(as.character(vdj_usage_df$count))
  
  circos.clear()
  tiff(paste("vdj_junction_gDNA_long_",id_specimen[i],".tiff",sep=""),res=300,h=1800,w=1800)
  circos.par(gap.degree = c(rep(2, length(unique(v))-1), 10,
                            rep(2, length(unique(d))-1), 10,
                            rep(2, length(unique(j))-1), 10))
  
  grid.col<-c(colorRampPalette(brewer.pal(6,"Reds"))(length(unique(v))),
              colorRampPalette(brewer.pal(6,"Greens"))(length(unique(d))),
              colorRampPalette(brewer.pal(6,"Purples"))(length(unique(j))))
  names(grid.col)<-c(as.character(unique(v)),as.character(unique(d)),as.character(unique(j)))
  
  chordDiagram(vdj_usage_df, grid.col = grid.col, annotationTrack = "grid",
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
}





###Circular plot for V and J junctionfor(i in 1:length(specimen)){
specimen<-unique(data_qc_gDNA_long_qc$specimen_label)
for(i in 1:length(specimen)){
  print(i)
  data_qc_gDNA_long_qc_specimen<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label==specimen[i]),]
  v_j_gene_mat<-unclass(table(data_qc_gDNA_long_qc_specimen$v_gene,data_qc_gDNA_long_qc_specimen$j_gene))
  v_j_gene_mat_sign<-v_j_gene_mat[na.omit(match(vgenes_sign,rownames(v_j_gene_mat))),]
  circos.clear()
  tiff(paste("v_j_junction_gDNA_long_",specimen[i],".tiff",sep=""),res=300,h=1800,w=1800)
  circos.par(gap.degree = c(rep(2, nrow(v_j_gene_mat_sign)-1), 10,
                            rep(2, ncol(v_j_gene_mat_sign)-1), 10))
  
  grid.col<-c(rev(colorRampPalette(brewer.pal(9,"Set1"))(nrow(v_j_gene_mat_sign))),
              colorRampPalette(brewer.pal(6,"Blues"))(ncol(v_j_gene_mat_sign)))
  
  names(grid.col)<-c(as.character(rownames(v_j_gene_mat_sign)),as.character(colnames(v_j_gene_mat_sign)))
  
  chordDiagram(v_j_gene_mat_sign, grid.col = grid.col, annotationTrack = "grid",
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
} 



###Using data frames

specimen<-unique(data_qc_gDNA_long_qc$specimen_label)
for(i in 1:length(specimen)){
  print(i)
  data_qc_gDNA_long_qc_specimen<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label==specimen[i]),]
  v_j_gene_df<-data.frame(table(data_qc_gDNA_long_qc_specimen$v_gene,data_qc_gDNA_long_qc_specimen$j_gene))
  v_j_gene_df.nonceros<-v_j_gene_df[which(v_j_gene_df$Freq!=0),]
  v_j_gene_df.nonceros_oder<-v_j_gene_df.nonceros[order(v_j_gene_df.nonceros[,1]), ]
  circos.clear()
  tiff(paste("v_j_junction_gDNA_long_",specimen[i],".tiff",sep=""),res=300,h=1800,w=1800)
  circos.par(gap.degree = c(rep(2, length(unique(v_j_gene_df.nonceros_oder[[1]]))-1), 10,
                            rep(2, length(unique(v_j_gene_df.nonceros_oder[[2]]))-1), 10))
  
  grid.col<-c(colorRampPalette(brewer.pal(9,"Set1"))(length(unique(v_j_gene_df.nonceros_oder[[1]]))),
              colorRampPalette(brewer.pal(6,"Blues"))(length(unique(v_j_gene_df.nonceros_oder[[2]]))))
  names(grid.col)<-c(as.character(unique(v_j_gene_df.nonceros_oder[[1]])),as.character(unique(v_j_gene_df.nonceros_oder[[2]])))
  
  chordDiagram(v_j_gene_df.nonceros_oder, grid.col = grid.col, annotationTrack = "grid",
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
} 

















####### J genes #########
jgenes<-as.data.frame(unclass(table(data_qc_gDNA_long_qc$specimen_label,data_qc_gDNA_long_qc$j_gene)))
id.spec<-match(rownames(jgenes),reads_clones_annot_Long_qc$specimen_id)
jgenes<-cbind(jgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec])
colnames(jgenes)[7:9]<-c("clin","time","Sample_id")

jusage<-matrix(NA,nrow(jgenes),(ncol(jgenes)-3))
for (i in 1:6){
  jusage[,i]<-jgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(jusage)<-colnames(jgenes)[1:6]
rownames(jusage)<-rownames(jgenes)

jgenes_filter<-jgenes[which(jgenes$time==0 | jgenes$time==6 | jgenes$time==24),]
jusage_filter<-jusage[which(jgenes$time==0 | jgenes$time==6 | jgenes$time==24),]

#### heatmap

annotation_col = data.frame(
  sample = factor(jgenes_filter$Sample_id),
  time = as.factor(jgenes_filter$time),
  clin = jgenes_filter$clin)

rownames(annotation_col)<-rownames(jusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
fill3=rainbow(26)
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),
                   sample= c("1005" = fill3[1],"1006" = fill3[2],"12001"= fill3[3],"12004"= fill3[4],"12008"= fill3[5],"38010"= fill3[6],
                             "38011"= fill3[7],"38012"= fill3[8],"38020"= fill3[9],"38021"=fill3[10], "41006"=fill3[11]  ,"89002"=fill3[12], "300003"=fill3[13] ,
                             "300005"=  fill3[14], "300006"= fill3[15],"300013"= fill3[16],"300015"= fill3[17],"301002"= fill3[18],"301003"= fill3[19],
                             "303001"= fill3[20], "303003"= fill3[21], "303005"= fill3[22], "303006"= fill3[23], "303007"= fill3[24], "303010"=fill3[25] , 
                             "303013"= fill3[26]))
tiff("heatmap_jgene_byclones.tiff",res=300,w=3000,h=2000)
pheatmap(t(jusage_filter),cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Purples"),border_color=F)
dev.off()


#########
## Statistical analysis to find significant jgenes
########
p<-NULL
for(i in 1:6){
   fm_null <- lmer(jusage[,i] ~ jgenes$clin + jgenes$time + (jgenes$time | jgenes$Sample_id),REML = F)
  fm_full <- lmer(jusage[,i] ~  jgenes$clin*jgenes$time + (jgenes$time | jgenes$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##non
  
} 

set.seed(2)
jusage_time<-data.frame(cbind(jusage,jgenes$time))
fit<-VSURF(x = jusage_time, y=factor(jgenes$clin),parallel = TRUE,ncores=4)
vars<-data.frame(jusage_time[,fit$varselect.interp])
set.seed(400)
rf_output <- randomForest(factor(jgenes$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=100)
tiff("MDSplot_jgene_byreads_RF.tiff",res=300,w=2000,h=2000)
MDSplot(rf_output, factor(jgenes$clin),palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "RF selection",cex=1.6)
legend("bottomleft", legend=levels(factor(vgenes$clin)),col=c("chartreuse4", "dodgerblue3","darkorange2"), pch = 20,pt.cex=1.6 ,cex=1.2)
dev.off()

jusage_sign<-vars
x<-t(scale(jusage_sign))
annotation_col = data.frame(
  clin = jgenes$clin)
rownames(annotation_col)<-rownames(jusage)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]))
tiff("heatmap_jgene_byclones.tiff",res=300,w=3000,h=2000)
pheatmap(t(jusage_sign),cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Purples"),border_color=F)
dev.off()


####### D genes #########
dgenes<-as.data.frame(unclass(table(data_qc_gDNA_long_qc$specimen_label,data_qc_gDNA_long_qc$d_gene)))
id.spec<-match(rownames(dgenes),reads_clones_annot_Long_qc$specimen_id)
dgenes<-cbind(dgenes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec])
colnames(dgenes)[32:34]<-c("clin","time","Sample_id")

dusage<-matrix(NA,nrow(dgenes),(ncol(dgenes)-3))
for (i in 1:31){
  dusage[,i]<-dgenes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(dusage)<-colnames(dgenes)[1:31]
rownames(dusage)<-rownames(dgenes)

dgenes_filter<-dgenes[which(dgenes$time==0 | dgenes$time==6 | dgenes$time==24),]
dusage_filter<-dusage[which(dgenes$time==0 | dgenes$time==6 | dgenes$time==24),]


#### heatmap
annotation_col = data.frame(
  sample = factor(dgenes_filter$Sample_id),
  time = as.factor(dgenes_filter$time),
  clin = dgenes_filter$clin)

rownames(annotation_col)<-rownames(dusage_filter)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
fill3=rainbow(26)
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),
                   sample= c("1005" = fill3[1],"1006" = fill3[2],"12001"= fill3[3],"12004"= fill3[4],"12008"= fill3[5],"38010"= fill3[6],
                             "38011"= fill3[7],"38012"= fill3[8],"38020"= fill3[9],"38021"=fill3[10], "41006"=fill3[11]  ,"89002"=fill3[12], "300003"=fill3[13] ,
                             "300005"=  fill3[14], "300006"= fill3[15],"300013"= fill3[16],"300015"= fill3[17],"301002"= fill3[18],"301003"= fill3[19],
                             "303001"= fill3[20], "303003"= fill3[21], "303005"= fill3[22], "303006"= fill3[23], "303007"= fill3[24], "303010"=fill3[25] , 
                             "303013"= fill3[26]))
tiff("heatmap_dgene_byclones.tiff",res=300,w=3000,h=2000)
pheatmap(t(dusage_filter),cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Greens"),border_color=F)
dev.off()



#########
## Statistical analysis to find significant dgenes
########
p<-NULL
for(i in 1:31){
  fm_null <- lmer(dusage[,i] ~ dgenes$clin + dgenes$time + (dgenes$time | dgenes$Sample_id),REML = F)
  fm_full <- lmer(dusage[,i] ~  dgenes$clin*dgenes$time + (dgenes$time | dgenes$Sample_id),REML = F)
  p[i]<-anova(fm_full, fm_null)[2,8] ##one
  
} 

set.seed(2)
dusage_time<-cbind(dusage,dgenes$time)
fit<-VSURF(x = dusage_time, y=factor(dgenes$clin),parallel = TRUE,ncores=4)
vars<-data.frame(dusage_time[,fit$varselect.interp])
set.seed(400)
rf_output <- randomForest(factor(dgenes$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=100)
tiff("MDSplot_dgene_byreads_RF.tiff",res=300,w=2000,h=2000)
MDSplot(rf_output, factor(dgenes$clin),palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "RF selection",cex=1.6)
legend("topleft", legend=levels(factor(dgenes$clin)),col=c("chartreuse4", "dodgerblue3","darkorange2"), pch = 20,pt.cex=1.6 ,cex=1.2)
dev.off()

dusage_sign<-vars
x<-t(scale(dusage_sign))
annotation_col = data.frame(
  clin = dgenes$clin)
rownames(annotation_col)<-rownames(dusage)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]))
tiff("heatmap_dgene_byreads_RF.tiff",res=300,w=3000,h=2000)
pheatmap(x,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"Greens"),border_color=F)
dev.off()

##### V J junction #######
data_qc_gDNA_long_qc$vj_genes<- paste(data_qc_gDNA_long_qc$v_gene,data_qc_gDNA_long_qc$j_gene,sep="_")
vj_genes<-as.data.frame(unclass(table(data_qc_gDNA_long_qc$specimen_label,data_qc_gDNA_long_qc$vj_genes)))
id.spec<-match(rownames(vj_genes),reads_clones_annot_Long_qc$specimen_id)
vj_genes<-cbind(vj_genes,reads_clones_annot_Long_qc$clin[id.spec],reads_clones_annot_Long_qc$time2[id.spec],reads_clones_annot_Long_qc$Sample_id[id.spec])
colnames(vj_genes)[326:328]<-c("clin","time","Sample_id")

vj_usage<-matrix(NA,nrow(vj_genes),(ncol(vj_genes)-3))
for (i in 1:325){
  vj_usage[,i]<-vj_genes[,i]/reads_clones_annot_Long_qc$clones_gDNA[id.spec]
}
colnames(vj_usage)<-colnames(vj_genes)[1:325]
rownames(vj_usage)<-rownames(vj_genes)

#### heatmap
x <- t(scale(vj_usage))

annotation_col = data.frame(
  clin = vj_genes$clin)
rownames(annotation_col)<-rownames(vj_usage)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]))
tiff("heatmap_vj_byclones.tiff",res=300,w=3000,h=2000)
pheatmap(x,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"RdBu"),border_color=F)
dev.off()

#########
## Statistical analysis to find significant jgenes
######## 


set.seed(2)
vjusage_time_clones<-cbind(vj_usage,vj_genes$time)
fit<-VSURF(x = vjusage_time_clones, y=factor(vj_genes$clin),parallel = TRUE,ncores=4)
vars<-data.frame(vjusage_time_clones[,fit$varselect.interp])
set.seed(400)
rf_output <- randomForest(factor(vj_genes$clin)~.,data=vars,proximity=TRUE, keep.forest=T,ntree=100)
tiff("MDSplot_vj_gene_byclones_RF.tiff",res=300,w=2000,h=2000)
MDSplot(rf_output, factor(jgenes$clin),palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "RF selection",cex=1.6)
legend("topleft", legend=levels(factor(vgenes$clin)),col=c("chartreuse4", "dodgerblue3","darkorange2"), pch = 20,pt.cex=1.6 ,cex=1.2)
dev.off()

vjusage_sign<-vars

x <- t(scale(vjusage_sign)) ##The scale is with the genes in the columns and samples in the rows

annotation_col = data.frame(
  clin = vj_genes$clin)
rownames(annotation_col)<-rownames(vj_usage)
fill=c("chartreuse4", "dodgerblue3","darkorange2")
ann_colors = list (clin = c(NP = fill[1], PNR = fill[2], PR = fill[3]))
tiff("heatmap_vjgene_byclones_RF.tiff",res=300,w=3000,h=2000)
pheatmap(x,cex=1.0,annotation_col = annotation_col,annotation_colors = ann_colors,rownames = F,color =  brewer.pal(6,"RdBu"),border_color=F)
dev.off()




