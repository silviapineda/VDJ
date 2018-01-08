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
setwd("/Users/Pinedasans/VDJ/SummaryResults/AllClones/")
load("~/VDJ/Data/clonesInferedAll_gDNA_90.Rdata")
clone_distribution<-colSums(clone_type_gDNA_num_reduced)
clone_distribution[order(clone_distribution,decreasing = T)]

###Filter out by clone distribution [all clones that only have 1 BCR in 1 sample]
clone_distribution_filter<-clone_distribution[which(clone_distribution>1)] 

barplot(clone_distribution_filter[order(clone_distribution_filter,decreasing = T)])

###Example clone with the maximum levels
table(data_gDNA_long_qc[grep("IGHV3-48_IGHJ4_30_3795",data_gDNA_long_qc$CloneId2),c("clin","time")])


####Filter the clone matrix to perform the analysis
id<-match(names(clone_distribution_filter),colnames(clone_type_gDNA_num_reduced))
clone_type_gDNA_num_filter<-clone_type_gDNA_num_reduced[,id]
clone_type_gDNA_num_filter2<-clone_type_gDNA_num_filter[,colSums(clone_type_gDNA_num_filter!=0)>3] ##194 clones

matrix_clones<-cbind(clone_type_gDNA_df[,2],clone_type_gDNA_num_filter2)
set.seed(112233)
fit<-VSURF(x = matrix_clones, y=factor(clone_type_gDNA_df[,1]),parallel = TRUE,ncores=4)
clones<-data.frame(matrix_clones[,fit$varselect.interp])

set.seed(1000)
rf_output <- randomForest(factor(clone_type_gDNA_df[,1])~.,data=clones,proximity=TRUE, keep.forest=T,ntree=1000)
tiff("MDSplot_cloneType_RF.tiff",res=300,w=2000,h=2000)
clone_type_gDNA_df$timeplot<-ifelse(clone_type_gDNA_df$time==0,1,
                                           ifelse(clone_type_gDNA_df$time==6,1.5,2))
MDSplot(rf_output, factor(clone_type_gDNA_df$clin),cex=clone_type_gDNA_df$timeplot,
        palette =c("chartreuse4", "dodgerblue3","darkorange2"),main = "clones")
legend("topleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

###Heatmap
clone_type_gDNA_df_filter<-clone_type_gDNA_df[which(clone_type_gDNA_df$time==0 |
                                                      clone_type_gDNA_df$time==6 |  
                                                      clone_type_gDNA_df$time==24),1:4]
clones<-clones[which(clone_type_gDNA_df$time==0 |
                       clone_type_gDNA_df$time==6 |  
                       clone_type_gDNA_df$time==24),]

splitop<-strsplit(as.character(clone_type_gDNA_df_filter$subject_id),"_")
subjects<-unlist(lapply(splitop, `[[`, 1))
numbers<-strsplit(subjects, "e")
individuals<-paste0("individual",unlist(lapply(numbers, `[[`, 2)))

annotation_col = data.frame(
  individual = factor(individuals),
  time = as.factor(clone_type_gDNA_df_filter$time),
  clin = clone_type_gDNA_df_filter$clin)


fill=c("chartreuse4", "dodgerblue3","darkorange2")
fill2=brewer.pal(3,"Set3")
n <- length(unique(individuals))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
fill3 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ann_colors = list (clin = c("NP" = fill[1], "PNR" = fill[2], "PR" = fill[3]),
                   time = c("0" = fill2[1],"6" = fill2[2], "24" = fill2[3]),  
                   individual= c("individual1" = fill3[1],"individual2" = fill3[2],"individual3"= fill3[3],"individual4"= fill3[4],"individual5"= fill3[5],"individual6"= fill3[6],
                                 "individual7"= fill3[7],"individual8"= fill3[8],"individual9"= fill3[9],"individual10"=fill3[10], "individual11"=fill3[11]  ,"individual12"=fill3[12], "individual13"=fill3[13] ,
                                 "individual14"=  fill3[14], "individual16"= fill3[15],"individual17"= fill3[16],"individual18"= fill3[17],"individual19"= fill3[18],"individual20"= fill3[19],
                                 "individual21"= fill3[20], "individual22"= fill3[21], "individual23"= fill3[22], "individual24"= fill3[23], "individual25"= fill3[24], "individual26"=fill3[25]))


colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_cloneType_RF.tiff",res=300,w=3500,h=2500)
pheatmap(t(clones),cex=1.0,annotation_col = annotation_col,
         annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

prefix<-"clone"
sufix<-seq(1,dim(clone_type_gDNA_filter2)[2])
my_names<-paste(prefix,sufix,sep=".")
colnames(clone_type_gDNA_filter2)<-my_names
paste(my_names,collapse ="+")

matrix_clones<-cbind(clone_type_gDNA_df_filter[,1:4],clone_type_gDNA_filter2)
paste(my_names,collapse ="+")

library(glmmLasso)
lm1 <- glmmLasso(clin ~ clone.1+clone.2+clone.3+clone.4+clone.5+clone.6+clone.7+clone.8+clone.9+clone.10+
                   clone.11+clone.12+clone.13+clone.14+clone.15+clone.16+clone.17+clone.18+clone.19+clone.20
                 +clone.21+clone.22+clone.23+clone.24+clone.25+clone.26+clone.27+clone.28+clone.29+clone.30+
                   clone.31+clone.32+clone.33+clone.34+clone.35+clone.36+clone.37+clone.38+clone.39+clone.40+
                   clone.41+clone.42+clone.43+clone.44+clone.45+clone.46+clone.47+clone.48+clone.49+clone.50+
                   clone.51+clone.52+clone.53+clone.54+clone.55+clone.56+clone.57+clone.58+clone.59+clone.60+
                   clone.61+clone.62+clone.63+clone.64+clone.65+clone.66+clone.67+clone.68+clone.69+clone.70+
                   clone.71+clone.72+clone.73, family = acat(),
                 rnd = list(time=~1+Sample_id),switch.NR=TRUE,control=list(print.iter=TRUE),
                 lambda=10,data=matrix_clones)
library(corrplot)
corrplot(cor(clone_type_gDNA_filter2))

p_value=NULL
for (i in 59:dim(matrix_clones)[2]){
  print(i)
  fm_null <- lmer(matrix_clones[,i] ~ matrix_clones$clin + matrix_clones$time + (1 | Sample_id),data=matrix_clones,REML = F)
  fm_full <- lmer(matrix_clones[,i]  ~  matrix_clones$clin*matrix_clones$time + (1 | Sample_id) ,data=matrix_clones,REML = F)
  p_value[i]<-anova(fm_full, fm_null)[2,8]
}

p <- ggplot(fm_full, aes(x = time, y =matrix_clones[,i], colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)

library(gamlr)
fit <- gamlr(matrix_clones[,5:ncol(matrix_clones)], matrix_clones$clin, gamma=10, lambda.min.ratio=0.1,
             standardize=FALSE, family="binomial")

