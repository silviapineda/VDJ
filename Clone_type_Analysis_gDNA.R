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
library(vcd)
library(gridGraphics)
library(gridExtra)


##################
### Analysis in the clone data matrix
#################
setwd("/Users/Pinedasans/VDJ/ResultsAllClones//")
load("~/VDJ/Data/clones_gDNA.Rdata")


###Clones individual analysis across outcome and time points
##Delete Individual 8
clone_type_gDNA_df_no8<-clone_type_gDNA_df[which(clone_type_gDNA_df$individual_id!="Individual8"),]
id_sample<-match(rownames(clone_type_gDNA_df_no8),rownames(clone_type_gDNA_num_reduced))
clone_type_gDNA_num_reduced_no8<-clone_type_gDNA_num_reduced[id_sample,]

#############
### Analysis 1: present/no-present
############
matrix_clones_presence<-apply(clone_type_gDNA_num_reduced_no8,1,function(x) ifelse(x==0,"No Present","Present"))
clone_type_gDNA_df_no8$time<-replace(clone_type_gDNA_df_no8$time,clone_type_gDNA_df_no8$time==12,6)
matrix_clones_presence_time0<-matrix_clones_presence[,which(clone_type_gDNA_df_no8$time==0)]
matrix_clones_presence_time6<-matrix_clones_presence[,which(clone_type_gDNA_df_no8$time==6)]
matrix_clones_presence_time24<-matrix_clones_presence[,which(clone_type_gDNA_df_no8$time==24)]

p_value_0=NULL
p_value_6=NULL
p_value_24=NULL
for(i in 1:dim(matrix_clones_presence_time0)[1]){
  print(i)
  tab<-table(matrix_clones_presence_time0[i,],clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==0)])
  if(dim(tab)[1]>1){
    p_value_0[i]=fisher.test(tab)$p.value
  }
  tab<-table(matrix_clones_presence_time6[i,],clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==6)])
  if(dim(tab)[1]>1){
    p_value_6[i]=fisher.test(tab)$p.value
  }
  tab<-table(matrix_clones_presence_time24[i,],clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==24)])
  if(dim(tab)[1]>1){
    p_value_24[i]=fisher.test(tab)$p.value
  }
}

##time0
matrix_clones_presence_significant_time0<-matrix_clones_presence_time0[which(p_value_0<0.05),] #8
results_time0<-list()
results_time0_cdr3<-list()
plots<-list()
for(i in 1:dim(matrix_clones_presence_significant_time0)[1]){
  tab<-table(matrix_clones_presence_significant_time0[i,],clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==0)])
  id<-match(names(which(matrix_clones_presence_significant_time0[i,]=="Present")),rownames(clone_type_gDNA_df_no8))
  id_clone<-match(rownames(matrix_clones_presence_significant_time0)[i],colnames(clone_type_gDNA_df_no8))
  id_clone2<-grep(rownames(matrix_clones_presence_significant_time0)[i],data_gDNA_long$V_J_lenghCDR3_CloneId)
  results_time0_cdr3[[i]]<-unique(as.character(data_gDNA_long[id_clone2,"cdr3_seq_aa_q"]))
  results_time0[[i]]<-clone_type_gDNA_df_no8[id,c(1:3,id_clone)]
  #clone_status<-matrix_clones_presence_significant_time0[i,]
  #clinical_outcome<-clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==0)]
  #tab_str<-structable(clone_status~clinical_outcome)
  #mosaic(tab_str,shade=T,main= rownames(matrix_clones_presence_significant_time0)[i], gp = shading_hcl, gp_args = list(interpolate = c(1, 1.8)))
  #plots[[i]]<-grid.grab()
  
}
grid.newpage()
tiff("mosaicplot_time0.tiff",res=300,h=5500,w=5500)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]], plots[[7]],
             plots[[8]],ncol=4)
dev.off()
p_value_0[which(p_value_0<0.05)]
names(results_time0)<-rownames(matrix_clones_presence_significant_time0)
cat(capture.output(print(results_time0), file="clones_fisher_time0.txt"))


##time6
matrix_clones_presence_significant_time6<-matrix_clones_presence_time6[which(p_value_6<0.05),] #4
results_time6<-list()
results_time6_cdr3<-list()
plots<-list()
for(i in 1:dim(matrix_clones_presence_significant_time6)[1]){
  tab<-table(matrix_clones_presence_significant_time6[i,],clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==6)])
  id<-match(names(which(matrix_clones_presence_significant_time6[i,]=="Present")),rownames(clone_type_gDNA_df_no8))
  id_clone<-match(rownames(matrix_clones_presence_significant_time6)[i],colnames(clone_type_gDNA_df_no8))
  id_clone2<-grep(rownames(matrix_clones_presence_significant_time6)[i],data_gDNA_long$V_J_lenghCDR3_CloneId)
  results_time6[[i]]<-clone_type_gDNA_df_no8[id,c(1:3,id_clone)]
  results_time6_cdr3[[i]]<-unique(as.character(data_gDNA_long[id_clone2,"cdr3_seq_aa_q"]))
  #clone_status<-matrix_clones_presence_significant_time6[i,]
  #clinical_outcome<-clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==6)]
  #tab_str<-structable(clone_status~clinical_outcome)
  #mosaic(tab_str,shade=T,main= rownames(matrix_clones_presence_significant_time6)[i], 
  #       gp = shading_hcl, gp_args = list(interpolate = c(1, 1.8)))
  #plots[[i]]<-grid.grab()
  
}
grid.newpage()
tiff("mosaicplot_time6.tiff",res=300,h=5500,w=5500)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],ncol=4)
dev.off()

p_value_6[which(p_value_6<0.05)]
names(results_time6)<-rownames(matrix_clones_presence_significant_time6)
cat(capture.output(print(results_time6), file="clones_fisher_time6.txt"))

##time24
matrix_clones_presence_significant_time24<-matrix_clones_presence_time24[which(p_value_24<0.05),] #21
results_time24<-list()
results_time24_cdr3<-list()
plots<-list()
for(i in 1:dim(matrix_clones_presence_significant_time24)[1]){
  tab<-table(matrix_clones_presence_significant_time24[i,],clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==24)])
  id<-match(names(which(matrix_clones_presence_significant_time24[i,]=="Present")),rownames(clone_type_gDNA_df_no8))
  id_clone<-match(rownames(matrix_clones_presence_significant_time24)[i],colnames(clone_type_gDNA_df_no8))
  id_clone2<-grep(rownames(matrix_clones_presence_significant_time24)[i],data_gDNA_long$V_J_lenghCDR3_CloneId)
  results_time24_cdr3[[i]]<-unique(as.character(data_gDNA_long[id_clone2,"cdr3_seq_aa_q"]))
  results_time24[[i]]<-clone_type_gDNA_df_no8[id,c(1:3,id_clone)]
  #clone_status<-matrix_clones_presence_significant_time24[i,]
  #clinical_outcome<-clone_type_gDNA_df_no8$clin[which(clone_type_gDNA_df_no8$time==24)]
  #tab_str<-structable(clone_status~clinical_outcome)
  #mosaic(tab_str,shade=T,main= rownames(matrix_clones_presence_significant_time24)[i],
  #         gp = shading_hcl, gp_args = list(interpolate = c(1, 1.8)))
  #plots[[i]]<-grid.grab()
 }

grid.newpage()
tiff("mosaicplot_time24.tiff",res=300,h=5500,w=5500)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]], plots[[7]],
             plots[[8]],plots[[9]],plots[[10]],plots[[11]],plots[[12]],plots[[13]],plots[[14]],
             plots[[15]],plots[[16]],plots[[17]],plots[[18]],plots[[19]],plots[[20]],plots[[21]],ncol=5)
dev.off()

p_value_24[which(p_value_24<0.05)]
names(results_time24)<-rownames(matrix_clones_presence_significant_time24)
cat(capture.output(print(results_time24), file="clones_fisher_time24.txt"))

#########
##Analysis 2: Measuring Clonal expansion
########
reads_clones_annot_qc<-reads_clones_annot[which(reads_clones_annot$clones_gDNA>100),]
id_sample<-match(rownames(clone_type_gDNA_num_reduced_no8),reads_clones_annot_qc$specimen_id)
total_clones<-reads_clones_annot_qc[na.omit(id_sample),"clones_gDNA"]

clone_type_gDNA_num_reduced_no8_qc<-clone_type_gDNA_num_reduced_no8[na.omit(match(reads_clones_annot_qc$specimen_id,rownames(clone_type_gDNA_num_reduced_no8))),]
clone_type_gDNA_df_no8_qc<-clone_type_gDNA_df_no8[na.omit(match(reads_clones_annot_qc$specimen_id,rownames(clone_type_gDNA_df_no8))),]
clone_matrix_freq<-apply(clone_type_gDNA_num_reduced_no8_qc,2,function(x) x/total_clones)
clone_matrix_freq_time0<-clone_matrix_freq[which(clone_type_gDNA_df_no8_qc$time==0),]
clone_matrix_freq_time6<-clone_matrix_freq[which(clone_type_gDNA_df_no8_qc$time==6),]
clone_matrix_freq_time24<-clone_matrix_freq[which(clone_type_gDNA_df_no8_qc$time==24),]

clone_matrix_freq_time0_NP<-clone_matrix_freq[which(clone_type_gDNA_df_no8_qc$time==0 & clone_type_gDNA_df_no8_qc$cline=="NP"),]
clone_matrix_freq_time0_PNR<-clone_matrix_freq[which(clone_type_gDNA_df_no8_qc$time==0 & clone_type_gDNA_df_no8_qc$cline=="PNR"),]
clone_matrix_freq_time0_PR<-clone_matrix_freq[which(clone_type_gDNA_df_no8_qc$time==0 & clone_type_gDNA_df_no8_qc$cline=="PR"),]

##Distribution of clones
clone_distribution_time0<-colSums(clone_matrix_freq_time0)
clone_distribution_time0[order(clone_distribution_time0,decreasing = T)]
barplot(clone_distribution_time0[order(clone_distribution_time0,decreasing = T)][1:100])

ggplot(data=clone_distribution_time0, aes(x=clone_distribution_time0, fill=clone_type_gDNA_df_no8_qc$clin)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c("chartreuse4", "dodgerblue3","darkorange2")) + 
  theme(legend.position = "none")

clone_distribution_time0_NP<-colSums(clone_matrix_freq_time0_NP)
clone_distribution_time0_PNR<-colSums(clone_matrix_freq_time0_PNR)
clone_distribution_time0_PR<-colSums(clone_matrix_freq_time0_PR)

plot.multi.dens( list(clone_distribution_time0_NP,clone_distribution_time0_PR))
plot(density(clone_distribution_time0_NP))

###Find if there is something longitudinally
###To measure which individuals have clones across data points
sample<-unique(clone_type_gDNA_df$individual_id)

#filter by the ones who has three time points to study the persistance
sample_filter<-sample[-c(5,6,7,15,17,18,21,27)]
persistance<-list()  
for (i in 1:length(sample_filter)){
  print(i)
  clone_matrix_sample<-clone_type_gDNA_num_reduced[which(clone_type_gDNA_df$individual_id==sample_filter[i]),]
  time0<-clone_matrix_sample[1,which(clone_matrix_sample[1,]!=0)]
  time6<-clone_matrix_sample[2,which(clone_matrix_sample[2,]!=0)]
  time24<-clone_matrix_sample[3,which(clone_matrix_sample[3,]!=0)]
  persistance[[i]]<-intersect(names(time0),names(time6))
  persistance[[i]]<-unique(c(persistance[[i]],intersect(names(time6),names(time24))))
  
}
names(persistance)<-sample_filter


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


###To obtain the clonal persistant 
clone_df<-data.frame(table(data_gDNA_long$V_J_lenghCDR3_CloneId,data_gDNA_long$specimen_label))
colnames(clone_df)<-c("clone","specimen","count")

clone_df$clin = reads_clones_annot[clone_df$specimen,1]
clone_df$time = reads_clones_annot[clone_df$specimen,3]
clone_df$individual = reads_clones_annot[clone_df$specimen,6]

clone_df$subject_id = reads_clones_annot[clone_df$specimen,5]
splitop<-strsplit(as.character(clone_df$subject_id),"_")
subject<-unlist(lapply(splitop, `[[`, 1))
clone_df$subject_id = subject

clone_df_noceros = clone_df[which(clone_df$count!=0),] #120,934

id<-NULL
persistance_clones<-unlist(persistance)
for(i in 1:length(persistance_clones)){
  print(i)
  id<-c(id,grep(persistance_clones[i],clone_df$clone))
}

clone_df_persistace<-clone_df[id,]
clone_df_persistace$time<-as.numeric(as.character(clone_df_persistace$time))

##Delete the outlier sample8 
clone_df_persistace<-clone_df_persistace[which(clone_df_persistace$individual!="Individual8"),]
##Time 32 into time 24
clone_df_persistace$time<-replace(clone_df_persistace$time,clone_df_persistace$time==32,24)

g1<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="NP"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24)) + scale_y_continuous(limits = c(0,300)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ individual) + labs(x = "time", y = "Clonal persistance")

g2<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="PNR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24)) + scale_y_continuous(limits = c(0,300)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ individual) + labs(x = "time", y = "Clonal persistance")

g3<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="PR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24)) + scale_y_continuous(limits = c(0,300)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ individual) + labs(x = "time", y = "Clonal persistance")

tiff("Clonal_persistant_gDNA_long.tiff",res=300,h=1700,w=2700)
multiplot(g1,g2,g3)
dev.off()

unique_clones<-unique(clone_df_persistace$clone)
persitance_clones<-clone_type_gDNA_num_filter[,na.omit(match(unique_clones,colnames(clone_type_gDNA_num_filter)))]


#
#################################################
###Find those clones that are expanded in PR ###
################################################
PR_24<-clone_type_gDNA_num_reduced[which(clone_type_gDNA_df[,1]=="PR" & clone_type_gDNA_df[,2]==24),]
id_specimen<-match(rownames(PR_24),rownames(clone_type_gDNA_num_reduced))
PR_clones_95<-list()
PR_specimen_nocero<-list()
j=1
for (i in id_specimen){
  PR_specimen<-clone_type_gDNA_num_reduced[i,]
  xx<-PR_specimen[which(PR_specimen!=0)]
  PR_specimen_nocero[[j]]<-xx
  #barplot(PR_specimen_nocero[order(PR_specimen_nocero,decreasing=T)])
  PR_specimen_nocero_95<-xx[which(xx>quantile(xx,c(0.95)))]
  PR_clones_95[[j]]<-PR_specimen_nocero_95
  j<-j+1
}
barplot(PR_clones_95[[3]][order(PR_clones_95[[3]],decreasing=T)])

###Common clones PR at time 24 
clones12<-intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[2]]))
clones13<-intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[3]]))
clones14<-intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[4]]))
clones15<-intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[5]]))
clones23<-intersect(names(PR_specimen_nocero[[2]]),names(PR_specimen_nocero[[3]]))
clones24<-intersect(names(PR_specimen_nocero[[2]]),names(PR_specimen_nocero[[4]]))
clones25<-intersect(names(PR_specimen_nocero[[2]]),names(PR_specimen_nocero[[5]]))
clones34<-intersect(names(PR_specimen_nocero[[3]]),names(PR_specimen_nocero[[4]]))
clones35<-intersect(names(PR_specimen_nocero[[3]]),names(PR_specimen_nocero[[5]]))
clones45<-intersect(names(PR_specimen_nocero[[4]]),names(PR_specimen_nocero[[5]]))

clone123<-intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[2]])),
                               names(PR_specimen_nocero[[3]])) #IGHV3-23_IGHJ4_30_5345
clone124<-intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[2]])),
                    names(PR_specimen_nocero[[4]]))#IGHV3-23_IGHJ4_45_17123
clone125<-intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[2]])),
                    names(PR_specimen_nocero[[5]]))
clone234<-intersect(intersect(names(PR_specimen_nocero[[4]]),names(PR_specimen_nocero[[2]])),
                    names(PR_specimen_nocero[[3]]))
clone235<-intersect(intersect(names(PR_specimen_nocero[[5]]),names(PR_specimen_nocero[[2]])),
                    names(PR_specimen_nocero[[3]]))#"IGHV3-33_IGHJ4_30_2804"
clone134<-intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[4]])),
                    names(PR_specimen_nocero[[3]]))
clone135<-intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[5]])),
                    names(PR_specimen_nocero[[3]]))
clone245<-intersect(intersect(names(PR_specimen_nocero[[2]]),names(PR_specimen_nocero[[5]])),
                    names(PR_specimen_nocero[[4]]))

##All zero
clone1234<-intersect(intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[2]])),
          names(PR_specimen_nocero[[3]])),names(PR_specimen_nocero[[4]]))
clone1235<-intersect(intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[2]])),
                               names(PR_specimen_nocero[[3]])),names(PR_specimen_nocero[[5]]))
clone2345<-intersect(intersect(intersect(names(PR_specimen_nocero[[5]]),names(PR_specimen_nocero[[2]])),
                               names(PR_specimen_nocero[[3]])),names(PR_specimen_nocero[[4]]))
clone1345<-intersect(intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[5]])),
                               names(PR_specimen_nocero[[3]])),names(PR_specimen_nocero[[4]]))
clone1245<-intersect(intersect(intersect(names(PR_specimen_nocero[[1]]),names(PR_specimen_nocero[[5]])),
                               names(PR_specimen_nocero[[2]])),names(PR_specimen_nocero[[4]]))

##Clones shared at time 24 for three individuals
data_gDNA_long_qc[grep("IGHV3-23_IGHJ4_30_5345",data_gDNA_long_qc$CloneId_CDR3),c("individual_id","clin","time","cdr3_seq_aa_q")]
data_gDNA_long_qc[grep("IGHV3-23_IGHJ4_45_17123",data_gDNA_long_qc$CloneId_CDR3),c("individual_id","clin","time","cdr3_seq_aa_q")]
data_gDNA_long_qc[grep("IGHV3-33_IGHJ4_30_2804",data_gDNA_long_qc$CloneId_CDR3),c("individual_id","clin","time","cdr3_seq_aa_q")]
