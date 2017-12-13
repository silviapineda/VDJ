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

setwd("/Users/Pinedasans/VDJ/SummaryResults/AllClones/")
data<-read.csv("/Users/Pinedasans/VDJ/Data/ClonesInferedAll.csv")

##Call clones from the dataset generate with all from python
clones<- unique(data[,c("specimen_label","CloneId","amplification_template")])
clones_mat<-data.matrix(table(clones$specimen_label,clones$amplification_template))

###Add the two new columns to the diversity file
diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity.csv",header=T)
id<-match(diversity$specimen_id,rownames(clones_mat))
diversity_all<-cbind(diversity,as.numeric(clones_mat[id,1]),as.numeric(clones_mat[id,2]))
colnames(diversity_all)[43:44]<-c("clones_all_cDNA","clones_all_gDNA")

##############
#### gDNA ####
#############
diversity_gDNA_qc<-diversity_all[which(diversity_all$reads_gDNA>=100),]
diversity_long_gDNA<-diversity_gDNA_qc[which(diversity_gDNA_qc$clin=="NP" | diversity_gDNA_qc$clin=="PNR" | diversity_gDNA_qc$clin=="PR"),]
diversity_long_gDNA$clin<-factor(diversity_long_gDNA$clin, levels=c("NP", "PNR", "PR"))


#####Plot both clones call
tiff("plot_all_clones_sample_clones_gDNA.tiff",h=1600,w=2000,res=300)
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[diversity_long_gDNA$clin]
plot(diversity_long_gDNA$clones_all_gDNA, diversity_long_gDNA$clones_gDNA,col = cols,pch=20,
     ylab = "Clones by Sample",xlab = "Clones by All")
legend("topleft",legend=c("NP","PNR","PR"), col=c("chartreuse4", "dodgerblue3","darkorange2"), 
       pch=20,cex=c(1.2))
dev.off()

###Diversity analysis
tiff("boxplot_All_clones_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==0),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]), 
                  diversity_long_gDNA$clones_all_gDNA[which(diversity_long_gDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==6),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]), 
                diversity_long_gDNA$clones_all_gDNA[which(diversity_long_gDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_gDNA[which(diversity_long_gDNA$time2==24),], 
            aes(factor(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]), 
                diversity_long_gDNA$clones_all_gDNA[which(diversity_long_gDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")
multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==0)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==6)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]))
summary(glm(diversity_long_gDNA$clones_gDNA[which(diversity_long_gDNA$time2==24)] ~ diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]))

###longitudinal data
fm_null <- lmer(diversity_long_gDNA$clones_all_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA,REML = F)
fm_full <- lmer(diversity_long_gDNA$clones_all_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_All_clones_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA$clones_all_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()

###############
####  cDNA  ###
##############
diversity_cDNA_qc<-diversity_all[which(diversity_all$clones_all_cDNA>=100),]
diversity_long_cDNA<-diversity_cDNA_qc[which(diversity_cDNA_qc$clin=="NP" | diversity_cDNA_qc$clin=="PNR" | diversity_cDNA_qc$clin=="PR"),]
diversity_long_cDNA$clin<-factor(diversity_long_cDNA$clin, levels=c("NP", "PNR", "PR"))

##plot all calls
tiff("plot_all_clones_sample_clones_cDNA.tiff",h=1600,w=2000,res=300)
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[diversity_long_cDNA$clin]
plot(diversity_long_cDNA$clones_all_cDNA, diversity_long_cDNA$clones_cDNA,col = cols,pch=20,
     ylab = "Clones by Sample",xlab = "Clones by All")
legend("topleft",legend=c("NP","PNR","PR"), col=c("chartreuse4", "dodgerblue3","darkorange2"), 
       pch=20,cex=c(1.2))
dev.off()

##Diversity analysis
#diversity_long_cDNA<-diversity_long_cDNA[which(diversity_long_cDNA$subject_id!="sample22_6"),]
tiff("boxplot_clones_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==0),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]), 
                diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==6),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]), 
                diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_cDNA[which(diversity_long_cDNA$time2==24),], 
            aes(factor(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]), 
                diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")
multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==0)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]))
summary(glm(diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==6)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]))
summary(glm(diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==24)] ~ diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]))

###Longitudinal data
fm_null <- lmer(diversity_long_cDNA$clones_all_cDNA ~ clin + time + (time | Sample_id),data=diversity_long_cDNA,REML = F)
fm_full <- lmer(diversity_long_cDNA$clones_all_cDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_cDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_cDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_cDNA$clones_all_cDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()


############################################
#### Study common clones across samples ###
##########################################
data_gDNA<-data[which(data$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]
id<-match(data_gDNA_long$specimen_label,diversity_long_gDNA$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

clone_type_matrix<-t(as.data.frame(unclass(table(data_gDNA_long_qc$CloneId,factor(data_gDNA_long_qc$specimen_label)))))
id.spec<-match(rownames(clone_type_matrix),diversity_long_gDNA$specimen_id)
clone_type_matrix<-cbind(as.character(diversity_long_gDNA$clin[id.spec]),diversity_long_gDNA$time2[id.spec],
                         diversity_long_gDNA$Sample_id[id.spec],as.character(diversity_long_gDNA$subject_id[id.spec]),clone_type_matrix)
colnames(clone_type_matrix)[1:4]<-c("clin","time","Sample_id","subject_id")

clone_type_matrix_num<-clone_type_matrix[,5:ncol(clone_type_matrix)]
clone_type_matrix_num2<-t(apply(clone_type_matrix_num,1,as.numeric))
colnames(clone_type_matrix_num2)<-colnames(clone_type_matrix[,5:ncol(clone_type_matrix)])

clone_type_matrix_reduced<-clone_type_matrix_num2[,which(colSums(clone_type_matrix_num2)!=0)]

clone_distribution<-colSums(clone_type_matrix_reduced)
clone_distribution[order(clone_distribution,decreasing = T)]




lm1 <- glmmLasso(clin ~ clone_type_matrix[,5:ncol(clone_type_matrix)], rnd = (time | Sample_id),
                 lambda=10, data = clone_type_matrix,family=acat())

