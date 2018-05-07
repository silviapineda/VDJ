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
## cDNA ##
##########
data_cDNA<-data_merge[which(data_merge$amplification_template=="cDNA"),]
data_cDNA_long<-data_cDNA[which(data_cDNA$clin=="NP" | data_cDNA$clin=="PNR" | data_cDNA$clin=="PR"),]


diversity<-read.csv("/Users/Pinedasans/VDJ/Data/diversity_AllClones_cDNA.csv",header=T)
reads_clones_annot_cDNA<-reads_clones_annot[which(reads_clones_annot$cDNA=="cDNA"),]
id<-match(reads_clones_annot_cDNA$specimen_id,diversity[,1])
diversity_reads_clones<-cbind(reads_clones_annot_cDNA,diversity[id,-1])

##Filter by those that has less than 1000 reads
diversity_cDNA<-diversity_reads_clones[which(diversity_reads_clones$reads_cDNA>1000),] #75 samples
min(diversity_cDNA$reads_cDNA) #62173
##Down-sampling to 62173 reads
id<-unlist(lapply(diversity_cDNA$specimen_id, function(x) grep(x,data_cDNA_long$specimen_label)))
data_cDNA_long_to_down<-data_cDNA_long[id,]
specimen_id<-unique(data_cDNA_long_to_down$specimen_label)
data_cDNA_long_downsampled<-NULL
for (i in specimen_id){
  print(i)
  data_cDNA_long_to_down_reads<-data_cDNA_long_to_down[which(data_cDNA_long_to_down$specimen_label==i),]
  id_down<-sample(1:dim(data_cDNA_long_to_down_reads)[1],62173)
  data_cDNA_long_downsampled<-rbind(data_cDNA_long_downsampled,data_cDNA_long_to_down_reads[id_down,])
}
save(data_cDNA_long_downsampled,file="/Users/Pinedasans/VDJ/Data/data_cDNA_long_downsampled.Rdata")

load("/Users/Pinedasans/VDJ/Data/data_cDNA_long_downsampled.Rdata")
##Prepare data for the analysis
data_cDNA_long_downsampled$specimen_label<-factor(data_cDNA_long_downsampled$specimen_label)
read_count <- table(data_cDNA_long_downsampled$specimen_label)
clones_count<- unique(data_cDNA_long_downsampled[,c("specimen_label","V_J_lenghCDR3_CloneId","isotype")])
clones<-data.matrix(table(clones_count$specimen_label,clones_count$isotype))

specimen_unique<-unique(data_cDNA_long_downsampled$specimen_label)
entropy_IGHA<-NULL
entropy_IGHD<-NULL
entropy_IGHG<-NULL
entropy_IGHM<-NULL
for (i in 1:length(specimen_unique)){
  print(i)
  data_specimen_unique<-data_cDNA_long_downsampled[which(data_cDNA_long_downsampled$specimen_label==specimen_unique[i]),]
  clones_specimen_IGHA<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHA"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHD<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHD"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHG<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHG"),"V_J_lenghCDR3_CloneId"]
  clones_specimen_IGHM<-data_specimen_unique[which(data_specimen_unique$isotype=="IGHM"),"V_J_lenghCDR3_CloneId"]
  
  fi_IGHA<-as.numeric(table(clones_specimen_IGHA))/length(clones_specimen_IGHA)
  fi_IGHD<-as.numeric(table(clones_specimen_IGHD))/length(clones_specimen_IGHD)
  fi_IGHG<-as.numeric(table(clones_specimen_IGHG))/length(clones_specimen_IGHG)
  fi_IGHM<-as.numeric(table(clones_specimen_IGHM))/length(clones_specimen_IGHM)
  
  hi_IGHA<-fi_IGHA*log2(fi_IGHA)
  hi_IGHD<-fi_IGHD*log2(fi_IGHD)
  hi_IGHG<-fi_IGHG*log2(fi_IGHG)
  hi_IGHM<-fi_IGHM*log2(fi_IGHM)
  
  entropy_IGHA[i]=-sum(hi_IGHA)
  entropy_IGHD[i]=-sum(hi_IGHD)
  entropy_IGHG[i]=-sum(hi_IGHG)
  entropy_IGHM[i]=-sum(hi_IGHM)
  
}


id<-match(specimen_unique,diversity_cDNA$specimen_id)

diversity_cDNA_downsampled<-data.frame(cbind(as.character(diversity_cDNA$specimen_id[id]),as.character(diversity_cDNA$clin[id]),
                                  as.numeric(diversity_cDNA$time[id]),as.character(diversity_cDNA$Individual.id[id]),read_count,clones[,c(2,3,5,6)],entropy_IGHA,
                                  entropy_IGHD,entropy_IGHG,entropy_IGHM))
colnames(diversity_cDNA_downsampled)<-c("specimen_id","clin","time","individual_id","reads","clones_IGHA","clones_IGHD",
                                        "clones_IGHG","clones_IGHM","entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")

diversity_cDNA_downsampled$time2<-replace(diversity_cDNA_downsampled$time,diversity_cDNA_downsampled$time==12,6)
diversity_cDNA_downsampled$clin<-factor(diversity_cDNA_downsampled$clin)
diversity_cDNA_downsampled$clones_IGHA<-as.numeric(as.character(diversity_cDNA_downsampled$clones_IGHA))
diversity_cDNA_downsampled$clones_IGHD<-as.numeric(as.character(diversity_cDNA_downsampled$clones_IGHD))
diversity_cDNA_downsampled$clones_IGHG<-as.numeric(as.character(diversity_cDNA_downsampled$clones_IGHG))
diversity_cDNA_downsampled$clones_IGHM<-as.numeric(as.character(diversity_cDNA_downsampled$clones_IGHM))
diversity_cDNA_downsampled$reads<-as.numeric(as.character(diversity_cDNA_downsampled$reads))
diversity_cDNA_downsampled$entropy_IGHA<-as.numeric(as.character(diversity_cDNA_downsampled$entropy_IGHA))
diversity_cDNA_downsampled$entropy_IGHD<-as.numeric(as.character(diversity_cDNA_downsampled$entropy_IGHD))
diversity_cDNA_downsampled$entropy_IGHG<-as.numeric(as.character(diversity_cDNA_downsampled$entropy_IGHG))
diversity_cDNA_downsampled$entropy_IGHM<-as.numeric(as.character(diversity_cDNA_downsampled$entropy_IGHM))
diversity_cDNA_downsampled$time<-as.numeric(as.character(diversity_cDNA_downsampled$time))

##################################
######Statistical Analysis ######
#################################
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
j=1
PNR_6_clones<-NULL
PNR_24_clones<-NULL
PR_6_clones<-NULL
PR_24_clones<-NULL
for (i in c("clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")){
  print(i)
  assign(paste0("diversity_cDNA_downsampled_",i),diversity_cDNA_downsampled[which(diversity_cDNA_downsampled[match(i,colnames(diversity_cDNA_downsampled))]>100),])
  tiff(paste0("boxplot_clones_cDNA_",i,".tiff"),h=1800,w=3000,res=300)
  p2 = ggplot(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),], 
              aes(factor(get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6)]), 
                  get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))],fill=clin)) + 
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clinical outcome", y = "Number of clones") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none") + ylim(0,20000) +  theme(text = element_text(size=15)) 
  
  p3 = ggplot(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),], 
              aes(factor(get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24)]), 
                  get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))],fill=clin)) +
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clinical outcome", y = "Number of clones") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none") + ylim(0,20000) + theme(text = element_text(size=15)) 
  
  grid.arrange(p2,p3,ncol=3)
  dev.off()
  PNR_6_clones[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                       get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6)])))[2,4]
  PR_6_clones[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                      get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6)])))[3,4]
  
  PNR_24_clones[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                       get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24)])))[2,4]
  PR_24_clones[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                      get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24)])))[3,4]
  
  j=j+1
}

###Entropy
j=1
PNR_6_entropy<-NULL
PNR_24_entropy<-NULL
PR_6_entropy<-NULL
PR_24_entropy<-NULL
clones<-c("clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
for (i in c("entropy_IGHA","entropy_IGHD","entropy_IGHG","entropy_IGHM")){
  print(i)
  assign(paste0("diversity_cDNA_downsampled_",i),diversity_cDNA_downsampled[which(diversity_cDNA_downsampled[match(clones[j],colnames(diversity_cDNA_downsampled))]>100),])
  tiff(paste0("boxplot_entropy_cDNA_downsampled",i,".tiff"),h=1800,w=3000,res=300)
  p2 = ggplot(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),], 
              aes(factor(get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6)]), 
                  get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))],fill=clin)) + 
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 6",x="Clinical outcome", y = "Shannon entropy") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none") + ylim(10,14) +  theme(text = element_text(size=15)) 
  
  p3 = ggplot(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),], 
              aes(factor(get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24)]), 
                  get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))],fill=clin)) +
    geom_violin(adjust = .5) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + labs(title="time 24",x="Clinical outcome", y = "Shannon entropy") + 
    stat_summary(fun.data=data_summary)  + theme(legend.position="none") + ylim(10,14) + theme(text = element_text(size=15)) 
  
  grid.arrange(p2,p3,ncol=3)
  dev.off()
  PNR_6_entropy[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                      get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6)])))[2,4]
  PR_6_entropy[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                     get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==6)])))[3,4]
  
  PNR_24_entropy[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                       get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24)])))[2,4]
  PR_24_entropy[j]<-coef(summary(glm(get(paste0("diversity_cDNA_downsampled_",i))[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24),match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))] ~ 
                                      get(paste0("diversity_cDNA_downsampled_",i))$clin[which(get(paste0("diversity_cDNA_downsampled_",i))$time2==24)])))[3,4]
  
  j=j+1
}


##Clones
clones<-c("clones_IGHA","clones_IGHD","clones_IGHG","clones_IGHM")
j<-1
p_value_clones<-NULL
for (i in clones){
  assign(paste0("diversity_cDNA_downsampled_",i),diversity_cDNA_downsampled[which(diversity_cDNA_downsampled[match(i,colnames(diversity_cDNA_downsampled))]>100),])
  fm_null <- lmer(get(paste0("diversity_cDNA_downsampled_",i))[match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))][,1] ~ 
                    clin + time + (1 | individual_id),data=get(paste0("diversity_cDNA_downsampled_",i)),REML = F)
  fm_full <- lmer(get(paste0("diversity_cDNA_downsampled_",i))[match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))][,1] ~ 
                    clin*time + (1 | individual_id),data=get(paste0("diversity_cDNA_downsampled_",i)),REML = F)
  
  p_value_clones[j]<-anova(fm_full, fm_null)[2,8] 
  j=j+1
  
  tiff(paste0("plot_lmer_clones_downsampled_cDNA_",i,".tiff"),h=1200,w=1400,res=300)
  p <- ggplot(fm_full, aes(x = time, y = get(paste0("diversity_cDNA_downsampled_",i))[match(i,colnames(get(paste0("diversity_cDNA_downsampled_",i))))][,1], 
                           colour = clin)) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
    geom_point(size=2) +
    geom_smooth(method="lm",size=1) +
    labs(x = "Time (months)",y = "Number of clones") + theme_bw() + theme_light()
  
  
  print(p)
  dev.off()
}
names(p_value_clones)<-clones


##entropy longitudinal
fm_null <- lmer(diversity_cDNA_downsampled$entropy ~ clin + time + (time | individual_id),data=diversity_cDNA_downsampled,REML = F)
fm_full <- lmer(diversity_cDNA_downsampled$entropy ~  clin*time + (time | individual_id) ,data=diversity_cDNA_downsampled,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_entropy_downsampled_cDNA.tiff",h=1200,w=1400,res=300)
p <- ggplot(fm_full, aes(x = time, y = diversity_cDNA_downsampled$entropy, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=1.2) +
  geom_smooth(method="lm",size=0.8) +
  labs(x = "Time (months)",y = "Number of clones") + theme_bw() + theme_light()

print(p)
dev.off()
