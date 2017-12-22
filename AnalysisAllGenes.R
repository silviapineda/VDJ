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

####Save the data processed
save(clones_mat,diversity_all,data,file="~/VDJ/Data/Data_all_90.Rdata")

load("~/VDJ/Data/Data_all_90.Rdata")

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

###Barplot with the clones
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
tiff("barplot_clones_all_gDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_gDNA$clones_all_gDNA[which(diversity_long_gDNA$time2==0)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==0)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==0)],
        cex.names=0.8,las=2,ylim = c(0,8000),ylab = c("Clones"))
legend(0, 4000, legend=levels(diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_gDNA$clones_all_gDNA[which(diversity_long_gDNA$time2==6)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==6)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==6)],
        cex.names=0.8,las=2,ylim = c(0,8000),ylab = c("Clones"))
barplot(diversity_long_gDNA$clones_all_gDNA[which(diversity_long_gDNA$time2==24)],
        col = COLOR[diversity_long_gDNA$clin[which(diversity_long_gDNA$time2==24)]],
        names.arg = diversity_long_gDNA$subject_id[which(diversity_long_gDNA$time2==24)],
        cex.names=0.8,las=2,ylim = c(0,8000),ylab = c("Clones"))
dev.off()

###Filter out the sample 8
diversity_long_gDNA_filter<-diversity_long_gDNA[which(diversity_long_gDNA$Sample_id!="301002"),]

###Diversity analysis
tiff("boxplot_All_clones_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_gDNA_filter[which(diversity_long_gDNA_filter$time2==0),], 
            aes(factor(diversity_long_gDNA_filter$clin[which(diversity_long_gDNA_filter$time2==0)]), 
                diversity_long_gDNA_filter$clones_all_gDNA[which(diversity_long_gDNA_filter$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_gDNA_filter[which(diversity_long_gDNA_filter$time2==6),], 
            aes(factor(diversity_long_gDNA_filter$clin[which(diversity_long_gDNA_filter$time2==6)]), 
                diversity_long_gDNA_filter$clones_all_gDNA[which(diversity_long_gDNA_filter$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_gDNA_filter[which(diversity_long_gDNA_filter$time2==24),], 
            aes(factor(diversity_long_gDNA_filter$clin[which(diversity_long_gDNA_filter$time2==24)]), 
                diversity_long_gDNA_filter$clones_all_gDNA[which(diversity_long_gDNA_filter$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")
multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_gDNA_filter$clones_gDNA[which(diversity_long_gDNA_filter$time2==0)] ~ diversity_long_gDNA_filter$clin[which(diversity_long_gDNA_filter$time2==0)]))
summary(glm(diversity_long_gDNA_filter$clones_gDNA[which(diversity_long_gDNA_filter$time2==6)] ~ diversity_long_gDNA_filter$clin[which(diversity_long_gDNA_filter$time2==6)]))
summary(glm(diversity_long_gDNA_filter$clones_gDNA[which(diversity_long_gDNA_filter$time2==24)] ~ diversity_long_gDNA_filter$clin[which(diversity_long_gDNA_filter$time2==24)]))

###longitudinal data
fm_null <- lmer(diversity_long_gDNA_filter$clones_all_gDNA ~ clin + time + (time | Sample_id),data=diversity_long_gDNA_filter,REML = F)
fm_full <- lmer(diversity_long_gDNA_filter$clones_all_gDNA ~  clin*time + (time | Sample_id) ,data=diversity_long_gDNA_filter,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_All_clones_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_gDNA_filter$clones_all_gDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()

###############
####  cDNA  ###
##############
diversity_cDNA_qc<-diversity_all[which(diversity_all$clones_all_cDNA>13941),]
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

##barplot
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
tiff("barplot_clones_all_cDNA.tiff",res=300,w=3000,h=2500)
par(mfrow=c(3,1))
barplot(diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==0)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==0)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==0)],
        cex.names=0.8,las=2,ylim = c(0,85000),ylab = c("Clones"))
legend(0.2, 82000, legend=levels(diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]),col=COLOR,pch=15, cex=1)
barplot(diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==6)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==6)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==6)],
        cex.names=0.8,las=2,ylim = c(0,85000),ylab = c("Clones"))
barplot(diversity_long_cDNA$clones_all_cDNA[which(diversity_long_cDNA$time2==24)],
        col = COLOR[diversity_long_cDNA$clin[which(diversity_long_cDNA$time2==24)]],
        names.arg = diversity_long_cDNA$subject_id[which(diversity_long_cDNA$time2==24)],
        cex.names=0.8,las=2,ylim = c(0,85000),ylab = c("Clones"))
dev.off()

##Diversity analysis

##Filter the sample 22
diversity_long_cDNA_filter<-diversity_long_cDNA[which(diversity_long_cDNA$subject_id!="sample22_6"),]

tiff("boxplot_clones_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(diversity_long_cDNA_filter[which(diversity_long_cDNA_filter$time2==0),], 
            aes(factor(diversity_long_cDNA_filter$clin[which(diversity_long_cDNA_filter$time2==0)]), 
                diversity_long_cDNA_filter$clones_all_cDNA[which(diversity_long_cDNA_filter$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Clones")
p2 = ggplot(diversity_long_cDNA_filter[which(diversity_long_cDNA_filter$time2==6),], 
            aes(factor(diversity_long_cDNA_filter$clin[which(diversity_long_cDNA_filter$time2==6)]), 
                diversity_long_cDNA_filter$clones_all_cDNA[which(diversity_long_cDNA_filter$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clin", y = "Clones")
p3 = ggplot(diversity_long_cDNA_filter[which(diversity_long_cDNA_filter$time2==24),], 
            aes(factor(diversity_long_cDNA_filter$clin[which(diversity_long_cDNA_filter$time2==24)]), 
                diversity_long_cDNA_filter$clones_all_cDNA[which(diversity_long_cDNA_filter$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clin", y = "Clones")
multiplot(p1,p2,p3)
dev.off()
summary(glm(diversity_long_cDNA_filter$clones_all_cDNA[which(diversity_long_cDNA_filter$time2==0)] ~ diversity_long_cDNA_filter$clin[which(diversity_long_cDNA_filter$time2==0)]))
summary(glm(diversity_long_cDNA_filter$clones_all_cDNA[which(diversity_long_cDNA_filter$time2==6)] ~ diversity_long_cDNA_filter$clin[which(diversity_long_cDNA_filter$time2==6)]))
summary(glm(diversity_long_cDNA_filter$clones_all_cDNA[which(diversity_long_cDNA_filter$time2==24)] ~ diversity_long_cDNA_filter$clin[which(diversity_long_cDNA_filter$time2==24)]))

###Longitudinal data
fm_null <- lmer(diversity_long_cDNA_filter$clones_all_cDNA ~ clin + time + (1 | Sample_id),data=diversity_long_cDNA_filter,REML = F)
fm_full <- lmer(diversity_long_cDNA_filter$clones_all_cDNA ~  clin*time + (1 | Sample_id) ,data=diversity_long_cDNA_filter,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_cDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_long_cDNA_filter$clones_all_cDNA, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()


############################################
#### Study common clones across gDNA    ###
##########################################
data_gDNA<-data[which(data$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]
id<-match(data_gDNA_long$specimen_label,diversity_long_gDNA_filter$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

##Build the matrix with the clones by samples
clone_type_gDNA<-t(as.data.frame(unclass(table(data_gDNA_long_qc$CloneId,factor(data_gDNA_long_qc$specimen_label)))))
id.spec<-match(rownames(clone_type_gDNA),diversity_long_gDNA$specimen_id)
clone_type_gDNA<-cbind(as.character(diversity_long_gDNA$clin[id.spec]),diversity_long_gDNA$time2[id.spec],
                       diversity_long_gDNA$Sample_id[id.spec],as.character(diversity_long_gDNA$subject_id[id.spec]),clone_type_gDNA)
colnames(clone_type_gDNA)[1:4]<-c("clin","time","Sample_id","subject_id")
clone_type_gDNA_df<-as.data.frame(clone_type_gDNA)

clone_type_gDNA_num<-clone_type_gDNA[,5:ncol(clone_type_gDNA)]
clone_type_gDNA_num2<-t(apply(clone_type_gDNA_num,1,as.numeric))
colnames(clone_type_gDNA_num2)<-colnames(clone_type_gDNA[,5:ncol(clone_type_gDNA)])

clone_type_gDNA_num_reduced<-clone_type_gDNA_num2[,which(colSums(clone_type_gDNA_num2)!=0)]
##147,617 clones that at least one sample has 
save(clone_type_gDNA_df,clone_type_gDNA_num_reduced,diversity_long_gDNA_filter,file="~/VDJ/Data/clonesInferedAll_gDNA_90.Rdata")


############################################
#### Study common clones across  cDNA ###
##########################################
data_cDNA<-data[which(data$amplification_template=="cDNA"),]
data_cDNA_long<-data_cDNA[which(data_cDNA$clin!="AR" & data_cDNA$clin!="pre-AR"),]
id<-match(data_cDNA_long$specimen_label,diversity_long_cDNA_filter$specimen_id)
data_cDNA_long_qc<-data_cDNA_long[which(is.na(id)==F),]
data_cDNA_long_qc$specimen_label<-factor(data_cDNA_long_qc$specimen_label)

##Build the matrix with the clones by samples
clone_type_cDNA<-t(as.data.frame(unclass(table(data_cDNA_long_qc$CloneId,factor(data_cDNA_long_qc$specimen_label)))))
id.spec<-match(rownames(clone_type_cDNA),diversity_long_cDNA$specimen_id)
clone_type_cDNA<-cbind(as.character(diversity_long_cDNA$clin[id.spec]),diversity_long_cDNA$time2[id.spec],
                       diversity_long_cDNA$Sample_id[id.spec],as.character(diversity_long_cDNA$subject_id[id.spec]),clone_type_cDNA)
colnames(clone_type_cDNA)[1:4]<-c("clin","time","Sample_id","subject_id")
clone_type_cDNA_df<-as.data.frame(clone_type_cDNA)

clone_type_cDNA_num<-clone_type_cDNA[,5:ncol(clone_type_cDNA)]
clone_type_cDNA_num2<-t(apply(clone_type_cDNA_num,1,as.numeric))
colnames(clone_type_cDNA_num2)<-colnames(clone_type_cDNA[,5:ncol(clone_type_cDNA)])

clone_type_cDNA_num_reduced<-clone_type_cDNA_num2[,which(colSums(clone_type_cDNA_num2)!=0)]
##2,355,583 clones that at least one sample has 
save(clone_type_cDNA_df,clone_type_cDNA_num_reduced,diversity_long_cDNA_filter,file="~/VDJ/Data/clonesInferedAll_cDNA_90.Rdata")



#####
## Common clones across gDNA anc cDNA ##
#####
load("~/VDJ/Data/clonesInferedAll_cDNA_90.Rdata")
load("~/VDJ/Data/clonesInferedAll_gDNA_90.Rdata")
specimen<-rownames(clone_type_gDNA_num_reduced)
common<-matrix(NA,length(specimen),length(specimen))
for (i in 1:length(specimen)){
  for (j in 1:length(specimen)){
    #gDNA
    clone_type_gDNA_reduced_specimen<-clone_type_gDNA_num_reduced[which(rownames(clone_type_gDNA_num_reduced)==specimen[i]),]
    clone_type_gDNA_reduced_specimen_present<-clone_type_gDNA_reduced_specimen[which(clone_type_gDNA_reduced_specimen!=0)]
    
    #cDNA
    clone_type_cDNA_reduced_specimen<-clone_type_cDNA_num_reduced[which(rownames(clone_type_cDNA_num_reduced)==specimen[j]),]
    clone_type_cDNA_reduced_specimen_present<-clone_type_cDNA_reduced_specimen[which(clone_type_cDNA_reduced_specimen!=0)]
    
    common[i,j]<-length(intersect(names(clone_type_gDNA_reduced_specimen_present),
                                  names(clone_type_cDNA_reduced_specimen_present)))
  }
}
write.csv(common,"common_90.csv",row.names = F)
common<-read.csv("common_90.csv")
colnames(common)<-specimen
rownames(common)<-specimen

common2<-common[which(colSums(common)!=0),]
common3<-common2[,which(colSums(common2)!=0)]

subject<-factor(diversity_long_gDNA_filter[match(colnames(common3),diversity_long_gDNA_filter$specimen_id),"subject_id"])
clones_gDNA<-diversity_long_gDNA_filter[match(colnames(common3),diversity_long_gDNA_filter$specimen_id),"clones_all_gDNA"]
clones_common_perc<-common3/clones_gDNA
colnames(clones_common_perc)<-subject
rownames(clones_common_perc)<-subject
#Plot the matrix
xval <- format(clones_common_perc, format="f", digits=1)
pal <- colorRampPalette(c(rgb(0.96,0.96,1), rgb(0.1,0.1,0.9)), space = "rgb")

library(gplots)
heatmap.2(clones_common_perc, Rowv=FALSE, Colv=FALSE, dendrogram="none", xlab="Columns", ylab="Rows", col=pal,cellnote=xval,
          tracecol="#303030", trace="none", notecol="black", 
          notecex=0.8, keysize = 1.5, margins=c(5, 5))

common_diag<-diag(as.matrix(common3))
common_diag_freq<-as.matrix(common_diag)/clones_gDNA

COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
names(common_diag_freq)<-subject
clin<-diversity_long_gDNA_filter[match(names(common_diag_freq),diversity_long_gDNA_filter$subject_id),"clin"]
barplot(common_diag_freq[order(common_diag_freq,decreasing = T)],col = COLOR[clin[order(common_diag_freq,decreasing = T)]], cex.names=0.8,las=2)

time<-diversity_long_gDNA_filter[match(names(common_diag_freq),diversity_long_gDNA_filter$subject_id),"time2"]

time_2<-time[which(time==6 | time==24)]
common_diag_freq_2<-common_diag_freq[which(time==6 | time==24)]
summary(glm(common_diag_freq_2~factor(time_2)))
library("RColorBrewer")
color<-brewer.pal(6,"Set1")
boxplot(common_diag_freq_2~factor(time_2),col=color[c(1,4)])

