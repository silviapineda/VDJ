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
load("~/VDJ/Data/clones_cDNA.Rdata")


###Filter by clones that at least are share in 2 samples
clone_type_cDNA_num_filter<-clone_type_cDNA_num_reduced[,colSums(clone_type_cDNA_num_reduced!=0)>2] #12,352
###present/no-present
matrix_clones_presence<-apply(clone_type_cDNA_num_filter,1,function(x) ifelse(x==0,"No Present","Present"))
clone_type_cDNA_df$time<-replace(clone_type_cDNA_df$time,clone_type_cDNA_df$time==12,6)
matrix_clones_presence_time6<-matrix_clones_presence[,which(clone_type_cDNA_df$time==6)]
matrix_clones_presence_time24<-matrix_clones_presence[,which(clone_type_cDNA_df$time==24)]

p_value_6=NULL
p_value_24=NULL
for(i in 1:dim(matrix_clones_presence_time6)[1]){
  print(i)
  tab<-table(matrix_clones_presence_time6[i,],clone_type_cDNA_df$clin[which(clone_type_cDNA_df$time==6)])
  if(dim(tab)[1]>1){
    p_value_6[i]=fisher.test(tab)$p.value
  }
  tab<-table(matrix_clones_presence_time24[i,],clone_type_cDNA_df$clin[which(clone_type_cDNA_df$time==24)])
  if(dim(tab)[1]>1){
    p_value_24[i]=fisher.test(tab)$p.value
  }
}

##time6
matrix_clones_presence_significant_time6<-matrix_clones_presence_time6[which(p_value_6<0.05),] #113
results_time6<-list()
plots<-list()
for(i in 1:dim(matrix_clones_presence_significant_time6)[1]){
  tab<-table(matrix_clones_presence_significant_time6[i,],clone_type_cDNA_df$clin[which(clone_type_cDNA_df$time==6)])
  id<-match(names(which(matrix_clones_presence_significant_time6[i,]=="Present")),rownames(clone_type_cDNA_df))
  id_clone<-match(rownames(matrix_clones_presence_significant_time6)[i],colnames(clone_type_cDNA_df))
  results_time6[[i]]<-clone_type_cDNA_df[id,c(1:3,id_clone)]
  # clone_status<-matrix_clones_presence_significant_time6[i,]
  # clinical_outcome<-clone_type_cDNA_df$clin[which(clone_type_cDNA_df$time==6)]
  # tab_str<-structable(clone_status~clinical_outcome)
  #mosaic(tab_str,shade=T,main= rownames(matrix_clones_presence_significant_time6)[i], 
         #gp = shading_hcl, gp_args = list(interpolate = c(1, 1.8)))
  #plots[[i]]<-grid.grab()
  
}
grid.newpage()
tiff("mosaicplot_time6_cDNA.tiff",res=300,h=5500,w=5500)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],ncol=4)
dev.off()

names(results_time6)<-rownames(matrix_clones_presence_significant_time6)
cat(capture.output(print(results_time6), file="clones_fisher_time6_cDNA.txt"))

results_time6_gDNA<-read.csv("clones_results_time6.csv")

##time24
matrix_clones_presence_significant_time24<-matrix_clones_presence_time24[which(p_value_24<0.05),] #21
results_time24<-list()
plots<-list()
for(i in 1:dim(matrix_clones_presence_significant_time24)[1]){
  tab<-table(matrix_clones_presence_significant_time24[i,],clone_type_cDNA_df$clin[which(clone_type_cDNA_df$time==24)])
  id<-match(names(which(matrix_clones_presence_significant_time24[i,]=="Present")),rownames(clone_type_cDNA_df))
  id_clone<-match(rownames(matrix_clones_presence_significant_time24)[i],colnames(clone_type_cDNA_df))
  results_time24[[i]]<-clone_type_cDNA_df[id,c(1:3,id_clone)]
  # clone_status<-matrix_clones_presence_significant_time24[i,]
  # clinical_outcome<-clone_type_cDNA_df$clin[which(clone_type_cDNA_df$time==24)]
  # tab_str<-structable(clone_status~clinical_outcome)
  # mosaic(tab_str,shade=T,main= rownames(matrix_clones_presence_significant_time24)[i],
  #        gp = shading_hcl, gp_args = list(interpolate = c(1, 1.8)))
  # plots[[i]]<-grid.grab()
}

grid.newpage()
tiff("mosaicplot_time24_cDNA.tiff",res=300,h=5500,w=5500)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]], plots[[7]],
             plots[[8]],plots[[9]],plots[[10]],plots[[11]],plots[[12]],plots[[13]],plots[[14]],
             plots[[15]],plots[[16]],plots[[17]],plots[[18]],plots[[19]],plots[[20]],plots[[21]],ncol=5)
dev.off()

names(results_time24)<-rownames(matrix_clones_presence_significant_time24)
cat(capture.output(print(results_time24), file="clones_fisher_time24_cDNA.txt"))
results_time24_gDNA<-read.csv("clones_results_time24.csv")

match(names(results_time6),results_time6_gDNA$X)

##Distribution of clones
clone_distribution<-colSums(clone_type_cDNA_num_filter)
clone_distribution[order(clone_distribution,decreasing = T)]

##Top 20 clones
barplot(clone_distribution[order(clone_distribution,decreasing = T)])

###Example clone with the maximum levels
table(data_cDNA_long[grep("IGHV4-59_IGHJ4_30_4700",data_cDNA_long$V_J_lenghCDR3_CloneId),c("clin","time")])
table(data_cDNA_long[grep("IGHV4-39_IGHJ5_51_803",data_cDNA_long$V_J_lenghCDR3_CloneId),c("specimen_label")])

matrix_clones<-cbind(clone_type_cDNA_df[,2],clone_type_cDNA_num_filter)
set.seed(112233)
fit<-VSURF(x = matrix_clones, y=factor(clone_type_cDNA_df[,1]),parallel = TRUE,ncores=4)
clones<-data.frame(matrix_clones[,fit$varselect.interp])

set.seed(1000)
rf_output <- randomForest(factor(clone_type_cDNA_df[,1])~.,data=clones,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_cloneType_RF_cDNA.tiff",res=300,w=2000,h=2000)
clone_type_cDNA_df$timeplot<-ifelse(clone_type_cDNA_df$time==0,1,
                                           ifelse(clone_type_cDNA_df$time==6,1.5,2))
MDSplot(rf_output, factor(clone_type_cDNA_df$clin),cex=clone_type_cDNA_df$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "clones")
legend("topleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

###Heatmap
clone_type_cDNA_df_filter<-clone_type_cDNA_df[which(clone_type_cDNA_df$time==0 |
                                                      clone_type_cDNA_df$time==6 |  
                                                      clone_type_cDNA_df$time==24),1:4]
clones<-clones[which(clone_type_cDNA_df$time==0 |
                       clone_type_cDNA_df$time==6 |  
                       clone_type_cDNA_df$time==24),]

splitop<-strsplit(as.character(clone_type_cDNA_df_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))
numbers<-strsplit(subjects, "e")
individuals<-paste0("individual",unlist(lapply(numbers, `[[`, 2)))

annotation_col = data.frame(
  individual = factor(individuals),
  time = factor(clone_type_cDNA_df_filter$time),
  clin = clone_type_cDNA_df_filter$clin)


fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- length(unique(individuals))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   individual= c("individual1" = fill3[1],"individual2" = fill3[2],"individual3"= fill3[3],"individual4"= fill3[4],"individual5"= fill3[5],"individual6"= fill3[6],
                                 "individual7"= fill3[7],"individual8"= fill3[8],"individual9"= fill3[9],"individual10"=fill3[10], "individual11"=fill3[11]  ,"individual12"=fill3[12], "individual13"=fill3[13] ,
                                 "individual14"=  fill3[14],"individual15"=  fill3[15], "individual16"= fill3[16],"individual17"= fill3[17],"individual18"= fill3[18],"individual19"= fill3[19],"individual20"= fill3[20],
                                 "individual21"= fill3[21], "individual22"= fill3[22], "individual23"= fill3[23], "individual24"= fill3[24], "individual25"= fill3[25], "individual26"=fill3[26],"individual27"=fill3[27]))

colfunc<-colorRampPalette(c("white","red4"))

tiff("heatmap_cloneType_RF_cDNA.tiff",res=300,w=3500,h=2500)
pheatmap(t(clones),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(1000),border_color=F)
dev.off()

write.csv(clones,file="clones_cDNA_RF.csv")
table(data_cDNA_long[grep(colnames(clones)[1],data_cDNA_long$V_J_lenghCDR3_CloneId),c("clin","time")])

##Study the clones selected
clones<-read.csv("/Users/Pinedasans/VDJ/ResultsAllClones/VgeneUsage-Clone/clones_cDNA_RF.csv")
id<-NULL
for(i in 1:31){
  print(i)
  id<-c(id,grep(colnames(clones)[i],data_cDNA_long$V_J_lenghCDR3_CloneId))
  
}
data_cDNA_long_clone<-data_cDNA_long[id,]
data_cDNA_long_clone$clin<-factor(data_cDNA_long_clone$clin)
data_cDNA_long_clone_barplot<-as.data.frame(table(data_cDNA_long_clone[,c("V_J_lenghCDR3_CloneId","sample_id","clin")]))

data_cDNA_long_clone_barplot$clin<-factor(data_cDNA_long_clone_barplot$clin)
data_cDNA_long_clone_barplot<-data_cDNA_long_clone_barplot[which(data_cDNA_long_clone_barplot$Freq!=0),]
data_cDNA_long_clone_barplot<-data_cDNA_long_clone_barplot[order(data_cDNA_long_clone_barplot$Freq),]

tiff("IndividualClonesRF.tiff",res=300,w=4000,h=2000)
ggplot(data_cDNA_long_clone_barplot[order(data_cDNA_long_clone_barplot$clin),], 
       aes(V_J_lenghCDR3_CloneId, y=Freq,fill=clin)) +
  geom_bar(stat="identity",position = "dodge",colour="black") + scale_fill_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_col(position = position_stack(reverse = TRUE),colour="black")
dev.off()

tiff("IndividualClonesRF.tiff",res=300,w=4000,h=2000)
ggplot(data_cDNA_long_clone_barplot[order(data_cDNA_long_clone_barplot$clin),], 
       aes(V_J_lenghCDR3_CloneId, y=Freq,fill=clin)) + 
  geom_bar(stat="identity",colour="black", position = "dodge",alpha=0.5) + 
  scale_fill_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

ggplot(df, aes(x=factor(names), y=values, fill=factor(num))) + 
  geom_bar(stat="identity", colour="black", position="dodge", alpha=0.1)



library(gamlr)
x<-data.matrix(clone_type_cDNA_num_filter)
y<-clone_type_cDNA_df[,1]
y<-factor(replace(y,y=="PNR","PR"))
cvfit <- cv.gamlr(x, y,family="binomial" , gamma=1, verb=TRUE)

fitlasso <- gamlr(x, y,family="binomial" , gamma=1, lambda.min.ratio=cvfit$lambda.min)

