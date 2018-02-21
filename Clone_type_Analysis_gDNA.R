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


###############
## Find if there is something longitudinally
##############

##Obtain the persistance by sample
##Each sample need to have at least two time points and we count the persistance across time points

###Delete the time 9 
clone_type_gDNA_df_no8_qc_notime9<-clone_type_gDNA_df_no8_qc[which(rownames(clone_type_gDNA_df_no8_qc)!="7_S63"),]
clone_type_gDNA_num_reduced_no8_qc_notime9<-clone_type_gDNA_num_reduced_no8_qc[which(rownames(clone_type_gDNA_num_reduced_no8_qc)!="7_S63"),]
###To measure which individuals have clones across data points
sample<-unique(clone_type_gDNA_df_no8_qc$individual_id)
#filter by the ones who has three time points to study the persistance
sample_filter<-sample[-c(5,7)]
persistance<-list()  
for (i in 1:length(sample_filter)){
  print(i)
  clone_matrix_sample<-clone_type_gDNA_num_reduced_no8_qc[which(clone_type_gDNA_df_no8_qc$individual_id==sample_filter[i]),]
  if(dim(clone_matrix_sample)[1]==3){
    time0<-clone_matrix_sample[1,which(clone_matrix_sample[1,]!=0)]
    time6<-clone_matrix_sample[2,which(clone_matrix_sample[2,]!=0)]
    time24<-clone_matrix_sample[3,which(clone_matrix_sample[3,]!=0)]
    persistance[[i]]<-intersect(names(time0),names(time6))
    persistance[[i]]<-unique(c(persistance[[i]],intersect(names(time6),names(time24))))
    persistance[[i]]<-unique(c(persistance[[i]],intersect(names(time0),names(time24))))
  } else {
    time0<-clone_matrix_sample[1,which(clone_matrix_sample[1,]!=0)]
    time6<-clone_matrix_sample[2,which(clone_matrix_sample[2,]!=0)]
    persistance[[i]]<-intersect(names(time0),names(time6))  
  }
}
names(persistance)<-sample_filter


###Number of clones that persist 
persistance_number<-NULL
for (i in 1:length(persistance)){
  persistance_number[i]<-length(persistance[[i]])
}
clinical_outcome<-c("NP","NP","NP","NP","NP","NP","NP","PNR","PNR","PNR","PNR","PNR","PNR","PNR","PNR","PNR",
                    "PR","PR","PR","PR","PR","PR","PR")


summary(glm(persistance_number~clinical_outcome))
boxplot(persistance_number~clinical_outcome)


##plotting persistence
clone_df<-data.frame(table(data_gDNA_long_qc$V_J_lenghCDR3_CloneId,data_gDNA_long_qc$specimen_label))
colnames(clone_df)<-c("clone","specimen","count")

colnames(reads_clones_annot)[4]<-c("specimen")
clone_merge<-merge(clone_df,reads_clones_annot[,c(1,3,4,6)],by = "specimen")

clone_df_noceros = clone_merge[which(clone_merge$count!=0),] #120,934

unique_clones<-unique(unlist(persistance))
id<-NULL
for (i in 1:length(unique_clones)){
  print(i)
  id<-c(id,grep(unique_clones[i],clone_df_noceros$clone))
}

clone_df_persistace<-clone_df_noceros[id,]
clone_df_persistace$time<-replace(clone_df_persistace$time,clone_df_persistace$time=="32","24")
clone_df_persistace$time<-as.numeric(as.character(clone_df_persistace$time))

g1<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="NP"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24)) + scale_y_continuous(limits = c(0,300)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ Individual.id) + labs(x = "time", y = "Clonal persistence")

g2<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="PNR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24)) + scale_y_continuous(limits = c(0,300)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ Individual.id) + labs(x = "time", y = "Clonal persistence")

g3<-ggplot(clone_df_persistace[which(clone_df_persistace$clin=="PR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24)) + scale_y_continuous(limits = c(0,300)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ Individual.id) + labs(x = "time", y = "Clonal persistence")

tiff("Clonal_persistant_gDNA_long.tiff",res=300,h=1700,w=2700)
multiplot(g1,g2,g3)
dev.off()


###Expanded clones
expansion_NP<-clone_df_persistace$count[which(clone_df_persistace$clin=="NP")]
expansion_PNR<-clone_df_persistace$count[which(clone_df_persistace$clin=="PNR")]
expansion_PR<-clone_df_persistace$count[which(clone_df_persistace$clin=="PR")]

expansion<-c(expansion_NP,expansion_PNR,expansion_PR)
clinical_outcome_expansion<-c(rep("NP",length(expansion_NP)),rep("PNR",length(expansion_PNR)),rep("PR",length(expansion_PR)))

summary(glm(expansion~clinical_outcome_expansion))
boxplot(expansion~clinical_outcome_expansion)

##Common clones that persist
library (plyr)
persistance_df <- ldply (persistance, data.frame)
colnames(persistance_df)<-c("Individual","clone")
persistance_df$clin<-c(rep("NP",55),rep("PNR",75),rep("PR",133))
###NP
NP_clones<-persistance_df[which(persistance_df$clin=="NP"),]
table(duplicated(NP_clones$clone))
##PNR
PNR_clones<-persistance_df[which(persistance_df$clin=="PNR"),]
table(duplicated(PNR_clones$clone))
individual_persistance_PNR<-list()
for(i in 1:9){
  individual_persistance_PNR[[i]]<-PNR_clones[grep(PNR_clones[duplicated(PNR_clones$clone),2][i],PNR_clones$clone),]
}
##PR
PR_clones<-persistance_df[which(persistance_df$clin=="PR"),]
table(duplicated(PR_clones$clone))
individual_persistance_PR<-list()
for(i in 1:8){
  individual_persistance_PR[[i]]<-PR_clones[grep(PR_clones[duplicated(PR_clones$clone),2][i],PR_clones$clone),]
}

results<-c(names(results_time0),names(results_time6),names(results_time24))
match(results,PR_clones[duplicated(PR_clones$clone),2])
PR_clones[duplicated(PR_clones$clone),2][c(3,7)]
#Match with previous resulsts for PR: IGHV3-23_IGHJ4_45_17123 IGHV4-34_IGHJ6_42_2263
match(results,PNR_clones[duplicated(PNR_clones$clone),2])
PNR_clones[duplicated(PNR_clones$clone),2][6]
##Match for PNR IGHV3-33_IGHJ4_30_2804