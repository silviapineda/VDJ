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

setwd("/Users/Pinedasans/VDJ/SummaryResults/ClonesByIndividual/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")


############################################
#### Study common clones across samples ###
##########################################

##gDNA

data_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]

diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity.csv",header=T)
diversity_gDNA_qc<-diversity[which(diversity$reads_gDNA>=100),]
diversity_long_gDNA<-diversity_gDNA_qc[which(diversity_gDNA_qc$clin=="NP" | diversity_gDNA_qc$clin=="PNR" | diversity_gDNA_qc$clin=="PR"),]
diversity_long_gDNA$clin<-factor(diversity_long_gDNA$clin, levels=c("NP", "PNR", "PR"))

id<-match(data_gDNA_long$specimen_label,diversity_long_gDNA$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

clone_type_matrix_gDNA<-t(as.data.frame(unclass(table(data_gDNA_long_qc$V_J_lenghCDR3_Clone_igh,factor(data_gDNA_long_qc$specimen_label)))))
id.spec<-match(rownames(clone_type_matrix_gDNA),diversity_long_gDNA$specimen_id)
clone_type_matrix_gDNA<-cbind(as.character(diversity_long_gDNA$clin[id.spec]),diversity_long_gDNA$time2[id.spec],
                              diversity_long_gDNA$Sample_id[id.spec],as.character(diversity_long_gDNA$subject_id[id.spec]),clone_type_matrix_gDNA)
colnames(clone_type_matrix_gDNA)[1:4]<-c("clin","time","Sample_id","subject_id")
clone_type_df_gDNA<-as.data.frame(clone_type_matrix_gDNA)

clone_type_matrix_num_gDNA<-clone_type_matrix_gDNA[,5:ncol(clone_type_matrix_gDNA)]
clone_type_matrix_num2_gDNA<-t(apply(clone_type_matrix_num_gDNA,1,as.numeric))
colnames(clone_type_matrix_num2_gDNA)<-colnames(clone_type_matrix_gDNA[,5:ncol(clone_type_matrix_gDNA)])

clone_type_matrix_reduced_gDNA<-clone_type_matrix_num2_gDNA[,which(colSums(clone_type_matrix_num2_gDNA)!=0)]
save(clone_type_df_gDNA,clone_type_matrix_reduced_gDNA,diversity_long_gDNA,file="~/VDJ/Data/clonesBySamples_gDNA.Rdata")


###To measure which individuals have clones across data points
sample<-unique(clone_type_df$Sample_id)
splitop<-strsplit(as.character(clone_type_df$subject_id),"_")
subject<-unique(unlist(lapply(splitop, `[[`, 1)))

#filter by the ones who has three time points to study the persistance
sample_filter<-sample[-c(3,5,6,7,11,16,17,20,26)]
subject_filter<-subject[-c(3,5,6,7,11,16,17,20,26)]
persistance<-list()  
for (i in 1:length(sample_filter)){
  print(i)
  clone_matrix_sample<-clone_type_matrix_reduced[which(clone_type_df$Sample_id==sample_filter[i]),]
  time0<-clone_matrix_sample[1,which(clone_matrix_sample[1,]!=0)]
  time6<-clone_matrix_sample[2,which(clone_matrix_sample[2,]!=0)]
  time24<-clone_matrix_sample[3,which(clone_matrix_sample[3,]!=0)]
  persistance[[i]]<-intersect(names(time0),names(time6))
  persistance[[i]]<-unique(c(persistance[[i]],intersect(names(time6),names(time24))))
  
}

######
##Fish plot for clones by time
#####
library(fishplot)
timepoints=c(0,6,24)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

for (i in 1:length(sample_filter)){
  print(i)
  if(length(persistance[[i]])>1){
    clone_matrix_sample<-clone_type_matrix_reduced[which(clone_type_df$Sample_id==sample_filter[i]),]
    clone_type_df_sample<-clone_type_df[which(clone_type_df$Sample_id==sample_filter[i]),1:4]
    id<-match(persistance[[i]],colnames(clone_matrix_sample))
    clone_fish<-clone_matrix_sample[,id]
    fish = createFishObject(t(clone_fish),parents = rep(0,ncol(clone_fish)),timepoints=timepoints,col=sample(col_vector, ncol(clone_fish)))
    fish = layoutClones(fish)
    tiff(paste0(clone_type_df_sample$clin[1],"_",subject_filter[i],".tiff"),res=300,h=1500,w=1500)
      fishPlot(fish,shape="spline",title.btm=paste0(clone_type_df_sample$clin[1],"_",subject_filter[i]),
           cex.title=0.5, vlines=c(0,6,24), 
           vlab=c("time 0","time 6","time 24"))
    dev.off()
  }
}

names(persistance)<-subject_filter
write.csv(persistance,file="persistance.csv")


###cDNA
data_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
data_cDNA_long<-data_cDNA[which(data_cDNA$clin!="AR" & data_cDNA$clin!="pre-AR"),]


diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity.csv",header=T)
diversity_cDNA_qc<-diversity[which(diversity$reads_cDNA>=100),]
diversity_long_cDNA<-diversity_cDNA_qc[which(diversity_cDNA_qc$clin=="NP" | diversity_cDNA_qc$clin=="PNR" | diversity_cDNA_qc$clin=="PR"),]
diversity_long_cDNA$clin<-factor(diversity_long_cDNA$clin, levels=c("NP", "PNR", "PR"))

id<-match(data_cDNA_long$specimen_label,diversity_long_cDNA$specimen_id)
data_cDNA_long_qc<-data_cDNA_long[which(is.na(id)==F),]
data_cDNA_long_qc$specimen_label<-factor(data_cDNA_long_qc$specimen_label)

clone_type_matrix_cDNA<-t(as.data.frame(unclass(table(data_cDNA_long_qc$V_J_lenghCDR3_Clone_igh,factor(data_cDNA_long_qc$specimen_label)))))
id.spec<-match(rownames(clone_type_matrix_cDNA),diversity_long_cDNA$specimen_id)
clone_type_matrix_cDNA<-cbind(as.character(diversity_long_cDNA$clin[id.spec]),diversity_long_cDNA$time2[id.spec],
                         diversity_long_cDNA$Sample_id[id.spec],as.character(diversity_long_cDNA$subject_id[id.spec]),clone_type_matrix_cDNA)
colnames(clone_type_matrix_cDNA)[1:4]<-c("clin","time","Sample_id","subject_id")
clone_type_df_cDNA<-as.data.frame(clone_type_matrix_cDNA)

clone_type_matrix_num_cDNA<-clone_type_matrix_cDNA[,5:ncol(clone_type_matrix_cDNA)]
clone_type_matrix_num2_cDNA<-t(apply(clone_type_matrix_num_cDNA,1,as.numeric))
colnames(clone_type_matrix_num2_cDNA)<-colnames(clone_type_matrix_cDNA[,5:ncol(clone_type_matrix_cDNA)])

clone_type_matrix_reduced_cDNA<-clone_type_matrix_num2_cDNA[,which(colSums(clone_type_matrix_num2_cDNA)!=0)]

save(clone_type_df_cDNA,clone_type_matrix_reduced_cDNA,diversity_long_cDNA,file="~/VDJ/Data/clonesBySamples_cDNA.Rdata")


#################################
### common clones by samples ###
###############################
load("~/VDJ/Data/clonesBySamples_gDNA.Rdata")
load("~/VDJ/Data/clonesBySamples_cDNA.Rdata")

specimen<-rownames(clone_type_matrix_reduced_gDNA)
common<-matrix(NA,length(specimen),length(specimen))
for (i in 1:length(specimen)){
  for (j in 1:length(specimen)){
    #gDNA
    clone_type_gDNA_reduced_specimen<-clone_type_matrix_reduced_gDNA[which(rownames(clone_type_matrix_reduced_gDNA)==specimen[i]),]
    clone_type_gDNA_reduced_specimen_present<-clone_type_gDNA_reduced_specimen[which(clone_type_gDNA_reduced_specimen!=0)]
    
    #cDNA
    clone_type_cDNA_reduced_specimen<-clone_type_matrix_reduced_cDNA[which(rownames(clone_type_matrix_reduced_cDNA)==specimen[j]),]
    clone_type_cDNA_reduced_specimen_present<-clone_type_cDNA_reduced_specimen[which(clone_type_cDNA_reduced_specimen!=0)]
    
    common[i,j]<-length(intersect(names(clone_type_gDNA_reduced_specimen_present),
                                  names(clone_type_cDNA_reduced_specimen_present)))
  }
}
write.csv(common,"common.csv")
common<-read.csv("common.csv")
colnames(common)<-specimen
rownames(common)<-specimen

common2<-common[which(colSums(common)!=0),]
common3<-common2[,which(colSums(common2)!=0)]

subject<-factor(diversity_long_gDNA[match(colnames(common3),diversity_long_gDNA$specimen_id),"subject_id"])
clones_gDNA<-diversity_long_gDNA[match(colnames(common3),diversity_long_gDNA$specimen_id),"clones_gDNA"]
clones_common_perc<-common3/clones_gDNA
colnames(clones_common_perc)<-subject
rownames(clones_common_perc)<-subject
#Plot the matrix
xval <- format(clones_common_perc, format="f", digits=2)
pal <- colorRampPalette(c(rgb(0.96,0.96,1), rgb(0.1,0.1,0.9)), space = "rgb")

library(gplots)
heatmap.2(as.matrix(clones_common_perc), Rowv=FALSE, Colv=FALSE, dendrogram="none", xlab="Columns", ylab="Rows", col=pal,cellnote=xval,
          tracecol="#303030", trace="none", notecol="black", 
          notecex=0.8, keysize = 1.5, margins=c(5, 5))

common_diag<-diag(as.matrix(common3))
common_diag_freq<-common_diag/clones_gDNA

COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
names(common_diag_freq)<-subject
clin<-diversity_long_gDNA[match(names(common_diag_freq),diversity_long_gDNA$subject_id),"clin"]
barplot(common_diag_freq[order(common_diag_freq,decreasing = T)],col = COLOR[clin[order(common_diag_freq,decreasing = T)]], cex.names=0.8,las=2)

time<-diversity_long_gDNA[match(names(common_diag_freq),diversity_long_gDNA$subject_id),"time2"]

time_2<-time[which(time==6 | time==24)]
common_diag_freq_2<-common_diag_freq[which(time==6 | time==24)]
summary(glm(common_diag_freq_2~factor(time_2)))
library("RColorBrewer")
color<-brewer.pal(6,"Set1")
boxplot(common_diag_freq_2~factor(time_2),col=color[c(1,4)])

