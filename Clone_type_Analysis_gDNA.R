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


##################
### Analysis in the clone data matrix
#################
setwd("/Users/Pinedasans/VDJ/ResultsAllClones//")
load("~/VDJ/Data/clones_gDNA.Rdata")

###Filter by clones that at least are share in 10% of the samples
clone_type_gDNA_num_filter<-clone_type_gDNA_num_reduced[,colSums(clone_type_gDNA_num_reduced!=0)>5] #103

clone_distribution<-colSums(clone_type_gDNA_num_reduced)
clone_distribution[order(clone_distribution,decreasing = T)]


##Top 20 clones
barplot(clone_distribution[order(clone_distribution,decreasing = T)][1:20])

###Example clone with the maximum levels
table(data_gDNA_long[grep("IGHV3-23_IGHJ6_45_3106",data_gDNA_long$V_J_lenghCDR3_CloneId),c("clin","time")])


###Find if there is something longitudinally
###To measure which individuals have clones across data points
sample<-unique(clone_type_gDNA_df$Sample_id)
splitop<-strsplit(as.character(clone_type_gDNA_df$subject_id),"_")
subject<-unique(unlist(lapply(splitop, `[[`, 1)))

#filter by the ones who has three time points to study the persistance
sample_filter<-sample[-c(5,6,7,15,17,18,21,27)]
subject_filter<-subject[-c(5,6,7,15,17,18,21,27)]
persistance<-list()  
for (i in 1:length(sample_filter)){
  print(i)
  clone_matrix_sample<-clone_type_gDNA_num_reduced[which(clone_type_gDNA_df$Sample_id==sample_filter[i]),]
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

for (i in 6:length(sample_filter)){
  print(i)
  if(length(persistance[[i]])>1){
    clone_matrix_sample<-clone_type_gDNA_num_reduced[which(clone_type_gDNA_df$Sample_id==sample_filter[i]),]
    clone_type_df_sample<-clone_type_gDNA_df[which(clone_type_gDNA_df$Sample_id==sample_filter[i]),1:4]
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

###To obtain the clonal persistant 
clone_df<-data.frame(table(data_gDNA_long$V_J_lenghCDR3_CloneId,data_gDNA_long$specimen_label))
colnames(clone_df)<-c("clone","specimen","count")

clone_df$clin = reads_clones_annot[clone_df$specimen,1]
clone_df$time = reads_clones_annot[clone_df$specimen,3]

clone_df$sample = reads_clones_annot[clone_df$specimen,2]

clone_df_noceros = clone_df[which(clone_df$count!=0),] #120,934

id<-NULL
persistance_clones<-unlist(persistance)
for(i in 1:length(persistance_clones)){
  print(i)
  id<-c(id,grep(persistance_clones[i],clone_df$clone))
}

clone_df_persistace<-clone_df[id,]
#clone_df_persistace$time<-as.numeric(as.character(clone_df_persistace$time))

g1<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="NP"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clonal persistance")

g2<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="PNR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clonal persistant")

g3<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="PR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clonal persistant")

tiff("Clonal_persistant_gDNA_long.tiff",res=300,h=1700,w=2700)
multiplot(g1,g2,g3)
dev.off()

####Individual analysis
matrix_clones<-cbind(clone_type_gDNA_df[,2],clone_type_gDNA_num_filter)
set.seed(112233)
fit<-VSURF(x = matrix_clones, y=factor(clone_type_gDNA_df[,1]),parallel = TRUE,ncores=4)
clones<-data.frame(matrix_clones[,fit$varselect.interp])

set.seed(1000)
rf_output <- randomForest(factor(clone_type_gDNA_df[,1])~.,data=clones,proximity=TRUE, keep.forest=T,ntree=1000)
###We do not see anything, but.. 
## Are the clones from cDNA overlapping with the ones in gDNA
clones_cDNA<-read.csv("clones_cDNA_RF.csv")
###None overlap with gDNA