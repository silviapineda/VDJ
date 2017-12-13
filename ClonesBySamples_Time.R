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
data_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]

diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity.csv",header=T)
diversity_gDNA_qc<-diversity[which(diversity$reads_gDNA>=100),]
diversity_long_gDNA<-diversity_gDNA_qc[which(diversity_gDNA_qc$clin=="NP" | diversity_gDNA_qc$clin=="PNR" | diversity_gDNA_qc$clin=="PR"),]
diversity_long_gDNA$clin<-factor(diversity_long_gDNA$clin, levels=c("NP", "PNR", "PR"))

id<-match(data_gDNA_long$specimen_label,diversity_long_gDNA$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

clone_type_matrix<-t(as.data.frame(unclass(table(data_gDNA_long_qc$V_J_lenghCDR3_Clone_igh,factor(data_gDNA_long_qc$specimen_label)))))
id.spec<-match(rownames(clone_type_matrix),diversity_long_gDNA$specimen_id)
clone_type_matrix<-cbind(as.character(diversity_long_gDNA$clin[id.spec]),diversity_long_gDNA$time2[id.spec],
                         diversity_long_gDNA$Sample_id[id.spec],as.character(diversity_long_gDNA$subject_id[id.spec]),clone_type_matrix)
colnames(clone_type_matrix)[1:4]<-c("clin","time","Sample_id","subject_id")
clone_type_df<-as.data.frame(clone_type_matrix)

clone_type_matrix_num<-clone_type_matrix[,5:ncol(clone_type_matrix)]
clone_type_matrix_num2<-t(apply(clone_type_matrix_num,1,as.numeric))
colnames(clone_type_matrix_num2)<-colnames(clone_type_matrix[,5:ncol(clone_type_matrix)])

clone_type_matrix_reduced<-clone_type_matrix_num2[,which(colSums(clone_type_matrix_num2)!=0)]


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
