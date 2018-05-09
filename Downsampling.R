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
library(ggplot2)
library(gridExtra)
library(grid)
library(lme4)

setwd("/Users/Pinedasans/VDJ/ResultsAllClones/")
load("/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")

##########
## gDNA ##
##########
data_gDNA<-data_merge[which(data_merge$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin=="NP" | data_gDNA$clin=="PNR" | data_gDNA$clin=="PR"),]


diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity_AllClones_gDNA.csv",header=T)
reads_clones_annot_gDNA<-reads_clones_annot[which(reads_clones_annot$gDNA=="gDNA"),]
id<-match(reads_clones_annot_gDNA$specimen_id,diversity[,1])
diversity_reads_clones<-cbind(reads_clones_annot_gDNA,diversity[id,-1])

##Filter by those that has less than 1000 reads
diversity_gDNA<-diversity_reads_clones[which(diversity_reads_clones$reads_gDNA>1000),] 
min(diversity_gDNA$reads_gDNA) #1062
##Down-sampling to 1062 reads
id<-unlist(lapply(diversity_gDNA$specimen_id, function(x) grep(x,data_gDNA_long$specimen_label)))
data_gDNA_long_to_down<-data_gDNA_long[id,]
specimen_id<-unique(data_gDNA_long_to_down$specimen_label)
data_gDNA_long_downsampled<-NULL
for (i in specimen_id){
  print(i)
  data_gDNA_long_to_down_reads<-data_gDNA_long_to_down[which(data_gDNA_long_to_down$specimen_label==i),]
  id_down<-sample(1:dim(data_gDNA_long_to_down_reads)[1],1062)
  data_gDNA_long_downsampled<-rbind(data_gDNA_long_downsampled,data_gDNA_long_to_down_reads[id_down,])
}
save(data_gDNA_long_downsampled,file="/Users/Pinedasans/VDJ/Data/data_gDNA_long_downsampled.Rdata")

##Prepare data for the analysis
data_gDNA_long_downsampled$specimen_label<-factor(data_gDNA_long_downsampled$specimen_label)
read_count <- table(data_gDNA_long_downsampled$specimen_label)
clones_count<- unique(data_gDNA_long_downsampled[,c("specimen_label","V_J_lenghCDR3_CloneId")])
clones<-data.matrix(table(clones_count$specimen_label))

specimen_unique<-unique(data_gDNA_long_downsampled$specimen_label)
entropy<-rep(NA,length(specimen_unique))
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_gDNA_long_downsampled[which(data_gDNA_long_downsampled$specimen_label==specimen_unique[i]),]
  clones_specimen<-data_specimen_unique[,"V_J_lenghCDR3_CloneId"]
  fi<-as.numeric(table(clones_specimen))/length(clones_specimen)
  hi<-fi*log2(fi)
  entropy[i]=-sum(hi) #entropy(table(clones_specimen)) returns the same result but by default is the natural logarithm (log)
}
names(entropy)<-specimen_unique

id<-match(specimen_unique,diversity_gDNA$specimen_id)

diversity_gDNA_downsampled<-data.frame(cbind(as.character(diversity_gDNA$specimen_id[id]),as.character(diversity_gDNA$clin[id]),
                                  as.numeric(diversity_gDNA$time[id]),as.character(diversity_gDNA$Individual.id[id]),read_count,clones[,1],entropy))
colnames(diversity_gDNA_downsampled)<-c("specimen_id","clin","time","individual_id","reads","clones","entropy")

diversity_gDNA_downsampled$time2<-replace(diversity_gDNA_downsampled$time,diversity_gDNA_downsampled$time==12,6)
diversity_gDNA_downsampled$clin<-factor(diversity_gDNA_downsampled$clin)
diversity_gDNA_downsampled$clones<-as.numeric(as.character(diversity_gDNA_downsampled$clones))
diversity_gDNA_downsampled$reads<-as.numeric(as.character(diversity_gDNA_downsampled$reads))
diversity_gDNA_downsampled$entropy<-as.numeric(as.character(diversity_gDNA_downsampled$entropy))
diversity_gDNA_downsampled$time<-as.numeric(as.character(diversity_gDNA_downsampled$time))
##################################
######Statistical Analysis ######
#################################
###Sample 8 is an outlier at 6 and 24
diversity_gDNA_downsampled<-diversity_gDNA_downsampled[which(diversity_gDNA_downsampled$individual_id!="Individual8"),]

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
tiff("boxplot_clones_gDNA_downsampled.tiff",h=1800,w=3000,res=300)
p1 = ggplot(diversity_gDNA_downsampled[which(diversity_gDNA_downsampled$time2==0),], 
            aes(factor(diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==0)]), 
                diversity_gDNA_downsampled$clones[which(diversity_gDNA_downsampled$time2==0)],fill=clin)) + ylim(0,1000) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +  theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Number of clones") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")

p2 = ggplot(diversity_gDNA_downsampled[which(diversity_gDNA_downsampled$time==6),], 
            aes(factor(diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time==6)]), 
                diversity_gDNA_downsampled$clones[which(diversity_gDNA_downsampled$time==6)],fill=clin)) + ylim(0,1000) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + theme(text = element_text(size=15)) +
  labs(title="time 6",x="Clinical outcome", y = "Number of clones")+ stat_summary(fun.data=data_summary) + theme(legend.position="none")


p3 = ggplot(diversity_gDNA_downsampled[which(diversity_gDNA_downsampled$time2==24),], 
            aes(factor(diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==24)]), 
                diversity_gDNA_downsampled$clones[which(diversity_gDNA_downsampled$time2==24)],fill=clin)) + ylim(0,1000) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +theme(text = element_text(size=15)) +
  labs(title="time 24",x="Clinical outcome", y = "Number of clones") + stat_summary(fun.data=data_summary)+ theme(legend.position="none")


#multiplot(p1,p2,p3)
grid.arrange(p1, p2,p3,ncol=3)
dev.off()

summary(glm(diversity_gDNA_downsampled$clones[which(diversity_gDNA_downsampled$time2==0)] ~ diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==0)]))
summary(glm(diversity_gDNA_downsampled$clones[which(diversity_gDNA_downsampled$time2==6)] ~ diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==6)]))
summary(glm(diversity_gDNA_downsampled$clones[which(diversity_gDNA_downsampled$time2==24)] ~ diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==24)]))

###Entropy
tiff("boxplot_entropy_gDNA_downsampled.tiff",h=1800,w=3000,res=300)
p1 = ggplot(diversity_gDNA_downsampled[which(diversity_gDNA_downsampled$time2==0),], 
            aes(factor(diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==0)]), 
                diversity_gDNA_downsampled$entropy[which(diversity_gDNA_downsampled$time2==0)],fill=clin)) + ylim(7,10) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +  theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Shannon entropy") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")

p2 = ggplot(diversity_gDNA_downsampled[which(diversity_gDNA_downsampled$time==6),], 
            aes(factor(diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time==6)]), 
                diversity_gDNA_downsampled$entropy[which(diversity_gDNA_downsampled$time==6)],fill=clin)) + ylim(7,10) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + theme(text = element_text(size=15)) +
  labs(title="time 6",x="Clinical outcome", y = "Shannon entropy")+ stat_summary(fun.data=data_summary) + theme(legend.position="none")


p3 = ggplot(diversity_gDNA_downsampled[which(diversity_gDNA_downsampled$time2==24),], 
            aes(factor(diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==24)]), 
                diversity_gDNA_downsampled$entropy[which(diversity_gDNA_downsampled$time2==24)],fill=clin)) + ylim(7,10) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +theme(text = element_text(size=15)) +
  labs(title="time 24",x="Clinical outcome", y = "Shannon entropy") + stat_summary(fun.data=data_summary)+ theme(legend.position="none")


#multiplot(p1,p2,p3)
grid.arrange(p1, p2,p3,ncol=3)
dev.off()

summary(glm(diversity_gDNA_downsampled$entropy[which(diversity_gDNA_downsampled$time2==0)] ~ diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==0)]))
summary(glm(diversity_gDNA_downsampled$entropy[which(diversity_gDNA_downsampled$time2==6)] ~ diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==6)]))
summary(glm(diversity_gDNA_downsampled$entropy[which(diversity_gDNA_downsampled$time2==24)] ~ diversity_gDNA_downsampled$clin[which(diversity_gDNA_downsampled$time2==24)]))



##Clones longitudinal
fm_null <- lmer(diversity_gDNA_downsampled$clones ~ clin + time + (time | individual_id),data=diversity_gDNA_downsampled,REML = F)
fm_full <- lmer(diversity_gDNA_downsampled$clones ~  clin*time + (time | individual_id) ,data=diversity_gDNA_downsampled,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clones_downsampled_gDNA.tiff",h=1200,w=1400,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_gDNA_downsampled$clones, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=1.2) +
  geom_smooth(method="lm",size=0.8) +
  labs(x = "Time (months)",y = "Number of clones") + theme_bw() + theme_light()

print(p)
dev.off()

##entropy longitudinal
fm_null <- lmer(diversity_gDNA_downsampled$entropy ~ clin + time + (time | individual_id),data=diversity_gDNA_downsampled,REML = F)
fm_full <- lmer(diversity_gDNA_downsampled$entropy ~  clin*time + (time | individual_id) ,data=diversity_gDNA_downsampled,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_entropy_downsampled_gDNA.tiff",h=1200,w=1400,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_gDNA_downsampled$entropy, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=1.2) +
  geom_smooth(method="lm",size=0.8) +
  labs(x = "Time (months)",y = "Number of clones") + theme_bw() + theme_light()

print(p)
dev.off()
