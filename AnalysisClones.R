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
library("ggplot2")
library(lme4)
library("RColorBrewer")
library("randomForest")
library("VSURF")
library(pheatmap)

setwd("/Users/Pinedasans/VDJ/ResultsAllClones/")
load("~/VDJ/Data/VDJ_clonesAllmerged.Rdata")

##Call clones from the dataset generate with all from python
clones<- unique(data_merge[,c("specimen_label","V_J_lenghCDR3_CloneId","amplification_template","isotype")])
clones$isotype2<-replace(clones$isotype,clones$amplification_template=="gDNA","gDNA")
clones$isotype2<-replace(clones$isotype2,clones$isotype2=="","UNMAPPED")
clones_mat<-data.matrix(table(clones$specimen_label,clones$isotype2))

############################################
#### Study common clones across gDNA    ###
##########################################
data_gDNA<-data_merge[which(data_merge$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin=="NP" | data_gDNA$clin=="PNR" | data_gDNA$clin=="PR"),]

##Build the matrix with the clones by samples
clone_type_gDNA<-t(as.data.frame(unclass(table(data_gDNA_long$V_J_lenghCDR3_CloneId,factor(data_gDNA_long$specimen_label)))))
id.spec<-match(rownames(clone_type_gDNA),reads_clones_annot$specimen_id)
clone_type_gDNA<-cbind(as.character(reads_clones_annot$clin[id.spec]),reads_clones_annot$time[id.spec],
                       as.character(reads_clones_annot$Individual.id[id.spec]),clone_type_gDNA)
colnames(clone_type_gDNA)[1:3]<-c("clin","time","individual_id")
clone_type_gDNA_df<-as.data.frame(clone_type_gDNA)

clone_type_gDNA_num<-clone_type_gDNA[,4:ncol(clone_type_gDNA)]
clone_type_gDNA_num2<-t(apply(clone_type_gDNA_num,1,as.numeric))
colnames(clone_type_gDNA_num2)<-colnames(clone_type_gDNA[,4:ncol(clone_type_gDNA)])

clone_type_gDNA_num_reduced<-clone_type_gDNA_num2[,which(colSums(clone_type_gDNA_num2)!=0)]
##118,223 clones that at least one sample has 




############################################
#### Study common clones across  cDNA ###
##########################################
data_cDNA<-data_merge[which(data_merge$amplification_template=="cDNA"),]
data_cDNA_long<-data_cDNA[which(data_cDNA$clin!="AR" & data_cDNA$clin!="pre-AR"),]
data_cDNA_long$isotype2<-replace(data_cDNA_long$isotype,data_cDNA_long$amplification_template=="gDNA","gDNA")
data_cDNA_long$isotype2<-replace(data_cDNA_long$isotype2,data_cDNA_long$isotype2=="","UNMAPPED")
data_cDNA_long$isotype2<-replace(data_cDNA_long$isotype2,data_cDNA_long$IGHM_naive_memory=="naive","IGHM_naive")
data_cDNA_long$isotype2<-replace(data_cDNA_long$isotype2,data_cDNA_long$IGHM_naive_memory=="memory","IGHM_memory")

##Build the matrix with the clones by samples
clone_type_cDNA<-t(as.data.frame(unclass(table(data_cDNA_long$V_J_lenghCDR3_CloneId,factor(data_cDNA_long$specimen_label)))))
id.spec<-match(rownames(clone_type_cDNA),reads_clones_annot$specimen_id)
clone_type_cDNA<-cbind(as.character(reads_clones_annot$clin[id.spec]),reads_clones_annot$time[id.spec],
                       as.character(reads_clones_annot$Individual.id[id.spec]),clone_type_cDNA)
colnames(clone_type_cDNA)[1:3]<-c("clin","time","individual_id")
clone_type_cDNA_df<-as.data.frame(clone_type_cDNA)

clone_type_cDNA_num<-clone_type_cDNA[,4:ncol(clone_type_cDNA)]
clone_type_cDNA_num2<-t(apply(clone_type_cDNA_num,1,as.numeric))
colnames(clone_type_cDNA_num2)<-colnames(clone_type_cDNA[,4:ncol(clone_type_cDNA)])

clone_type_cDNA_num_reduced<-clone_type_cDNA_num2[,which(colSums(clone_type_cDNA_num2)!=0)]
##2,419,988 clones that at least one sample has 
save(clone_type_cDNA_df,clone_type_cDNA_num_reduced,reads_clones_annot,data_cDNA_long,file="~/VDJ/Data/clones_cDNA.Rdata")


######
###Check for correlation with clones by indiviudal definition
######
clones_igh<- unique(data_merge[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template")])
clones_igh<-data.matrix(table(clones_igh$specimen_label,clones_igh$amplification_template))
colnames(clones_igh)<-c("clones_igh_cDNA","clones_igh_gDNA")

reads_clones_annot_long<-reads_clones_annot[which(reads_clones_annot$clin!="pre-AR" & reads_clones_annot$clin!="AR"),]
id_common<-match(reads_clones_annot_long$specimen_id,rownames(clones_igh))
clones_igh<-data.matrix(clones_igh[id_common,])

#####Plot both clones call
tiff("plot_clonesAll_clonesIgh_gDNA.tiff",h=1600,w=2000,res=300)
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[factor(reads_clones_annot_long$clin)]
plot(reads_clones_annot_long$clones_gDNA, clones_igh[,2],col = cols,pch=20,
     ylab = "Clones by Sample",xlab = "Clones by All")
legend("topleft",legend=c("NP","PNR","PR"), col=c("chartreuse4", "dodgerblue3","darkorange2"), 
       pch=20,cex=c(1.2))
dev.off()

tiff("plot_clonesAll_clonesIgh_cDNA.tiff",h=1600,w=2000,res=300)
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[factor(reads_clones_annot_long$clin)]
plot(reads_clones_annot_long$clones_cDNA, clones_igh[,1],col = cols,pch=20,
     ylab = "Clones by Sample",xlab = "Clones by All")
legend("topleft",legend=c("NP","PNR","PR"), col=c("chartreuse4", "dodgerblue3","darkorange2"), 
       pch=20,cex=c(1.2))
dev.off()

########################################
## Common clones across gDNA and cDNA ##
#######################################
load("~/VDJ/Data/clones_gDNA.Rdata")
load("~/VDJ/Data/clones_cDNA.Rdata")
specimen_common<-intersect(rownames(clone_type_gDNA_num_reduced),rownames(clone_type_cDNA_num_reduced))
common<-matrix(NA,length(specimen_common),length(specimen_common))
for (i in 1:length(specimen_common)){
  print(i)
  for (j in 1:length(specimen_common)){
    #gDNA
    clone_type_gDNA_reduced_specimen<-clone_type_gDNA_num_reduced[which(rownames(clone_type_gDNA_num_reduced)==specimen_common[i]),]
    clone_type_gDNA_reduced_specimen_present<-clone_type_gDNA_reduced_specimen[which(clone_type_gDNA_reduced_specimen!=0)]
    
    #cDNA
    clone_type_cDNA_reduced_specimen<-clone_type_cDNA_num_reduced[which(rownames(clone_type_cDNA_num_reduced)==specimen_common[j]),]
    clone_type_cDNA_reduced_specimen_present<-clone_type_cDNA_reduced_specimen[which(clone_type_cDNA_reduced_specimen!=0)]
    
    common[i,j]<-length(intersect(names(clone_type_gDNA_reduced_specimen_present),
                                  names(clone_type_cDNA_reduced_specimen_present)))
  }
}
write.csv(common,"common.csv",row.names = F)
common<-data.matrix(read.csv("common.csv"))

colnames(common)<-specimen_common
rownames(common)<-specimen_common

#QC
reads_clones_annot_qc<-reads_clones_annot[which(reads_clones_annot$clones_gDNA>100),]
id<-match(reads_clones_annot_qc$specimen_id,colnames(common))
common_qc<-common[na.omit(id),na.omit(id)]

subject<-factor(reads_clones_annot_qc[match(colnames(common_qc),reads_clones_annot_qc$specimen_id),"subject_id"])
clones_gDNA<-reads_clones_annot_qc[match(colnames(common_qc),reads_clones_annot_qc$specimen_id),"clones_gDNA"]

clones_common_perc<-common_qc/clones_gDNA
colnames(clones_common_perc)<-subject
rownames(clones_common_perc)<-subject
#Plot the matrix
xval <- format(clones_common_perc, format="f", digits=1)
pal <- colorRampPalette(c(rgb(0.96,0.96,1), rgb(0.1,0.1,0.9)), space = "rgb")

library(gplots)
heatmap.2(clones_common_perc, Rowv=FALSE, Colv=FALSE, dendrogram="none", xlab="Columns", ylab="Rows", col=pal,cellnote=xval,
          tracecol="#303030", trace="none", notecol="black", 
          notecex=0.8, keysize = 1.5, margins=c(5, 5))

common_diag<-diag(as.matrix(common_qc))
common_diag_freq<-as.matrix(common_diag)/clones_gDNA
common_diag_freq<-common_diag_freq[,1]

COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
names(common_diag_freq)<-subject
clin<-factor(reads_clones_annot[match(names(common_diag_freq),reads_clones_annot$subject_id),"clin"])
barplot(common_diag_freq[order(common_diag_freq,decreasing = T)],col = COLOR[clin[order(common_diag_freq,decreasing = T)]], 
        cex.names=0.8,las=2)

time<-reads_clones_annot[match(names(common_diag_freq),reads_clones_annot$subject_id),"time"]

time_2<-time[which(time==6 | time==24)]
common_diag_freq_2<-common_diag_freq[which(time==6 | time==24)]
summary(glm(common_diag_freq_2~factor(time_2)))
library("RColorBrewer")
color<-brewer.pal(6,"Set1")
boxplot(common_diag_freq_2~factor(time_2),col=color[c(1,4)])

