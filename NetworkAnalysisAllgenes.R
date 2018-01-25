rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies
###
### CITATION: 
###
### PROCESS: Network Analysis
###           
### DESCRIP: Analysis of the VDJ data
###         
###
### Author: Silvia Pineda
### Date: June, 2017
############################################################################################
library(ggplot2)
library("igraph")
library("ineq")
library("dplyr")
library(lme4)

setwd("/Users/Pinedasans/VDJ/ResultsAllClones//")
load("~/VDJ/Data/VDJ_clonesAllmerged.Rdata")


############
### gDNA ###
############
data_gDNA<-data[which(data$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]
id<-match(data_gDNA_long$specimen_label,diversity_long_gDNA$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

###########################################
###### To study the persistances clones ###
##The vertex are the BCRs and the frequency of each one
##The edges shows how many are forming a clone
##Example: Number of BCRs in 7_S01 is 795 (vertex) and number of edges is 6. 795-6 = 789 number of clones

##First: obtain the number of vertex
specimen<-unique(data_qc_gDNA_long_qc$specimen_label)
for (i in specimen){
  print(i)
  data_qc_gDNA_long_qc_specimen<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label==i),]
  df_specimen<-data_qc_gDNA_long_qc_specimen[,c("clone_CDR3","igh_clone_id")]
  groups <- group_by(df_specimen,clone_CDR3,igh_clone_id)
  #edges<-unique(data.frame(groups))
  assign(paste("vertex",i,sep=""),data.frame(table(data_qc_gDNA_long_qc_specimen$clone_CDR3)))
}

##Second: Obtain if there are BCR's (vertex) share across time
sample<-unique(data_qc_gDNA_long_qc$sample_id)
#filter by the ones who has three time points to study the persistance
sample_filter<-sample[-c(3,5,6,7,11,16,17,26)]
persistance<-list()
j=1
for (i in sample_filter){
  print(i)
  specimen<-reads_clones_annot_Long_qc$specimen_id[which(reads_clones_annot_Long_qc$Sample_id==i)]
  #persistance[[j]]<-intersect(intersect(get(paste("vertex",specimen[1],sep=""))[,1],get(paste("vertex",specimen[2],sep=""))[,1]),get(paste("vertex",specimen[3],sep=""))[,1]) #No persistance clones across 3 time points
  #persistance[[j]]<-intersect(get(paste("vertex",specimen[1],sep=""))[,1],get(paste("vertex",specimen[2],sep=""))[,1]) #Persistance from time 0 to time 24
  persistance[[j]]<-intersect(get(paste("vertex",specimen[2],sep=""))[,1],get(paste("vertex",specimen[3],sep=""))[,1]) #Persistance from time 6 to 24
  j=j+1
}
sample_subject_ID<-read.csv("/Users/Pinedasans/VDJ/Data/Sample_subject_ID.csv")
subject<-sample_subject_ID[match(sample_filter,sample_subject_ID$Sample_id),]
names(persistance)<-subject$subject_id




###Obtain vertex and edges to apply the nucleotides-assembly-1.0.jar made by Mikel
specimen<-unique(data_qc_gDNA_long_qc$specimen_label)
for (i in specimen){
  print(i)
  data_qc_gDNA_long_qc_specimen<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label==i),]
  df_specimen<-data_qc_gDNA_long_qc_specimen[,c("clone_CDR3","igh_clone_id")]
  groups <- group_by(df_specimen,clone_CDR3,igh_clone_id)
  edges<-unique(data.frame(groups))
  vertex<-data.frame(table(data_qc_gDNA_long_qc_specimen$clone_CDR3))
  write.table(edges,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/long_gDNA/edges",i,".txt",sep=""),sep="\t",row.names = F)
  write.table(vertex,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/vertex",i,".txt",sep=""),sep="\t",row.names = F)
}


##Obtain the vertex and cluster Gini index
vertex_max<-NULL
vertex_gini<-NULL
cluster_max<-NULL
cluster_gini<-NULL
num_reads_max_cluster<-NULL
clusters<-NULL
j<-1
specimen<-unique(data_qc_gDNA_long_qc$specimen_label)
for (i in specimen){
  assign(paste("edges",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/long_gDNA/edges",i,".txt",sep=""),header=T))
  assign(paste("vertex",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/vertex",i,".txt",sep=""),header=T))
  vertex_max[j]<-max(get(paste("vertex",i,sep=""))$Freq)
  vertex_gini[j]<-Gini(get(paste("vertex",i,sep=""))$Freq)
  cluster_max[j]<-max(table(get(paste("edges",i,sep=""))$igh_clone_id))
  clusters[j]<-sum(table(table(get(paste("edges",i,sep=""))$igh_clone_id)))
  num_reads_max_cluster[j]<-tail(table(table(get(paste("edges",i,sep=""))$igh_clone_id)),1)
  cluster_gini[j]<-Gini(table(get(paste("edges",i,sep=""))$igh_clone_id))
  j=j+1
}

id<-match(specimen, reads_clones_annot_Long_qc$specimen_id)
reads_clones_annot_Long_qc$cluster_gini[id]<-cluster_gini
reads_clones_annot_Long_qc$vertex_gini[id]<-vertex_gini
reads_clones_annot_Long_qc$vertex_max[id]<-vertex_max
reads_clones_annot_Long_qc$cluster_max[id]<-cluster_max
reads_clones_annot_Long_qc$num_reads_max_cluster[id]<-num_reads_max_cluster
reads_clones_annot_Long_qc$clusters[id]<-clusters
reads_clones_annot_Long_qc$clonal_expansion<-(reads_clones_annot_Long_qc$num_reads_max_cluster/reads_clones_annot_Long_qc$reads_gDNA)*100

reads_clones_annot_Long_qc$clin<-factor(reads_clones_annot_Long_qc$clin)

###Delete sample 8 form the analysis
reads_clones_annot_Long_qc<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$subject_id.1!="sample8_6" & reads_clones_annot_Long_qc$subject_id.1!="sample8_24"),]
                                                                             

tiff("plot_vertex_gini_cluster_gini_gDNA.tiff",h=1600,w=2000,res=300)
reads_clones_annot_Long_qc$timeplot<-ifelse(reads_clones_annot_Long_qc$time2==0,1,
                               ifelse(reads_clones_annot_Long_qc$time2==6,1.5,2))
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[reads_clones_annot_Long_qc$clin]
plot(cluster_gini, vertex_gini,cex=reads_clones_annot_Long_qc$timeplot,
        col = cols,pch=20,ylab = "Vextex Gini",xlab = "Cluster Gini")
legend("topleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()

tiff("plot_vertex_gini_cluster_gini_gDNA_1by1.tiff",h=2000,w=1800,res=300)
reads_clones_annot_Long_qc_time0<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==0),]
p1<-qplot(cluster_gini, vertex_gini, data= reads_clones_annot_Long_qc_time0, colour = clin) +
  scale_colour_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) +
  scale_x_continuous(limit = c(0,0.1)) + scale_y_continuous(limit = c(0,0.6)) + labs(title="time 0")

reads_clones_annot_Long_qc_time6<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==6),]
p2<-qplot(cluster_gini, vertex_gini, data= reads_clones_annot_Long_qc_time6, colour = clin) +
  scale_colour_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) +
  scale_x_continuous(limit = c(0,0.1)) + scale_y_continuous(limit = c(0,0.6)) + labs(title="time 6")

reads_clones_annot_Long_qc_time24<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==24),]
p3<-qplot(cluster_gini, vertex_gini, data= reads_clones_annot_Long_qc_time24, colour = clin) +
  scale_colour_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) +
  scale_x_continuous(limit = c(0,0.1)) + scale_y_continuous(limit = c(0,0.6)) + labs(title="time 24")
multiplot(p1,p2,p3)
dev.off()

anova(lm(cluster_gini~reads_clones_annot_Long_qc$clin))
tiff("boxplot_vertexGini_clin_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(reads_clones_annot_Long_qc$clusters~reads_clones_annot_Long_qc$clin,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="Degree of clonal expansion")
dev.off()

anova(lm(vertex_gini~reads_clones_annot_Long_qc$clin))
tiff("boxplot_vertexGini_clin_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(reads_clones_annot_Long_qc$vertex_gini~reads_clones_annot_Long_qc$clin,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="Degree of clonal expansion")
dev.off()


tiff("boxplot_vertex_gini_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==0),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time2==0)]), 
                reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Vertex Gini")
p2 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==6),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time2==6)]), 
                reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 6",x="Clin", y = "Vertex Gini")
p3 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==24),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time2==24)]), 
                reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 24",x="Clin", y = "Vertex Gini")
multiplot(p1,p2,p3)
dev.off()

summary(glm(reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==0)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==0)]))
summary(glm(reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==6)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==6)]))
summary(glm(reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==24)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==24)]))


tiff("boxplot_cluster_gini_gDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==0),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time2==0)]), 
                reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time2==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Cluster Gini")
p2 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==6),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time2==6)]), 
                reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time2==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 6",x="Clin", y = "Cluster Gini")
p3 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==24),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time2==24)]), 
                reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time2==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 24",x="Clin", y = "Cluster Gini")
multiplot(p1,p2,p3)
dev.off()

summary(glm(reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==0)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==0)]))
summary(glm(reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==6)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==6)]))
summary(glm(reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==24)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==24)]))

##Gini
fm_null <- lmer(reads_clones_annot_Long_qc$cluster_gini ~ clin + time + (time | Sample_id),data=reads_clones_annot_Long_qc,REML = F)
fm_full <- lmer(reads_clones_annot_Long_qc$cluster_gini ~  clin*time + (time | Sample_id) ,data=reads_clones_annot_Long_qc,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clusterGini_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = reads_clones_annot_Long_qc$cluster_gini, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()

fm_null <- lmer(reads_clones_annot_Long_qc$vertex_gini ~ clin + time + (time | Sample_id),data=reads_clones_annot_Long_qc,REML = F)
fm_full <- lmer(reads_clones_annot_Long_qc$vertex_gini ~  clin*time + (time | Sample_id) ,data=reads_clones_annot_Long_qc,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_vertexGini_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = reads_clones_annot_Long_qc$vertex_gini, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()

###Once is in network format I plot the network
##NP
specimen_NP<- c("7_S01", "7_S02", "7_S03", "7_S04", "7_S05", "7_S06", "7_S07", "7_S09", "7_S10", "7_S11", "7_S12", 
                "7_S14", "7_S17", "7_S18", "7_S21", "7_S22", "7_S23", "7_S24", "7_S25", "7_S26", "7_S27",
                "7_S28", "7_S29", "7_S30")
specimen_PNR<-c("7_S32", "7_S33", "7_S34", "7_S35", "7_S36", "7_S37", "7_S38", "7_S39", "7_S40", "7_S41", "7_S42",
                "7_S46", "7_S47", "7_S48" ,"7_S50" , "7_S51", "7_S53", "7_S54", "7_S55", "7_S56", "7_S57","7_S58","7_S59", "7_S60")

specimen_PR<-c("7_S61", "7_S62", "7_S63", "7_S64", "7_S65", "7_S66", "7_S67", "7_S68", "7_S69", "7_S70", "7_S71", "7_S72",
               "7_S73", "7_S74", "7_S75", "7_S76", "7_S77", "7_S78", "7_S79", "7_S82", "7_S83")

for(i in specimen_NP) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("chartreuse4")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}
specimen_PNR<-c("7_S58","7_S59", "7_S60")
for(i in specimen_PNR) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("dodgerblue3")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}


for(i in specimen_PR) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("darkorange2")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}

tiff("/Users/Pinedasans/VDJ/Results/Network/NP.tiff",res=300, width = 3000, height =6000)
par(mfrow=c(6,3))
files <- list.files("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/NP/")
for(i in files) {
  cat(i, "\n")
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/NP/",i, sep = ""))
  j<-substr(i,6,10)
  nodes<-get(paste("vertex",j,sep=""))
  net<-graph_from_data_frame(d=edges,vertices = nodes,directed=T)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("chartreuse4")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  
  #tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_",j,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  #dev.off()
}
dev.off()

#PNR
tiff("/Users/Pinedasans/VDJ/Results/Network/PNR.tiff",res=300, width = 3000, height =6000)
par(mfrow=c(6,3))
files <- list.files("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/PNR/")
for(i in files) {
  cat(i, "\n")
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/PNR/",i, sep = ""))
  j<-substr(i,6,10)
  nodes<-get(paste("vertex",j,sep=""))
  net<-graph_from_data_frame(d=edges,vertices = nodes,directed=T)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("dodgerblue3")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  #tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_",j,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  #dev.off()
}
dev.off()

##PR
tiff("/Users/Pinedasans/VDJ/Results/Network/PR.tiff",res=300, width = 3000, height =6000)
par(mfrow=c(6,3))
files <- list.files("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/PR/")
for(i in files) {
  cat(i, "\n")
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_gDNA/PR/",i, sep = ""))
  j<-substr(i,6,10)
  nodes<-get(paste("vertex",j,sep=""))
  net<-graph_from_data_frame(d=edges,vertices = nodes,directed=T)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("darkorange2")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  #tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_",j,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  #dev.off()
}
dev


##############
### AR #######
#############
data_qc_gDNA_AR<-data_qc_gDNA[which(data_qc_gDNA$clin=="AR" | data_qc_gDNA$clin=="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
reads_clones_annot_AR_qc<-reads_clones_annot_AR[which(reads_clones_annot_AR$reads_gDNA>100),]

id<-match(data_qc_gDNA_AR$specimen_label,reads_clones_annot_AR_qc$specimen_id)
data_qc_gDNA_AR_qc<-data_qc_gDNA_AR[which(is.na(id)==F),]
data_qc_gDNA_AR_qc$specimen_label<-factor(data_qc_gDNA_AR_qc$specimen_label)

data_qc_gDNA_AR_qc$clone_CDR3<-paste(data_qc_gDNA_AR_qc$igh_clone_id, data_qc_gDNA_AR_qc$cdr3_seq_aa_q,sep="")


###Save to apply the nucleotides-assembly-1.0.jar made by Mikel
specimen<-unique(data_qc_gDNA_AR_qc$specimen_label)
for (i in specimen){
  print(i)
  data_qc_gDNA_AR_qc_specimen<-data_qc_gDNA_AR_qc[which(data_qc_gDNA_AR_qc$specimen_label==i),]
  df_specimen<-data_qc_gDNA_AR_qc_specimen[,c("clone_CDR3","igh_clone_id")]
  groups <- group_by(df_specimen,clone_CDR3,igh_clone_id)
  edges<-unique(data.frame(groups))
  vertex<-data.frame(table(data_qc_gDNA_AR_qc_specimen$clone_CDR3))
  write.table(edges,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/AR_gDNA/edges",i,".txt",sep=""),sep="\t",row.names = F)
  write.table(vertex,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_gDNA//vertex",i,".txt",sep=""),sep="\t",row.names = F)
}

##Obtain the vertex and cluster Gini index
vertex_max<-NULL
vertex_gini<-NULL
cluster_max<-NULL
cluster_gini<-NULL
num_reads_max_cluster<-NULL
clusters<-NULL
j<-1
specimen<-unique(data_qc_gDNA_AR_qc$specimen_label)
for (i in specimen){
  assign(paste("edges",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/AR_gDNA/edges",i,".txt",sep=""),header=T))
  assign(paste("vertex",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_gDNA/vertex",i,".txt",sep=""),header=T))
  vertex_max[j]<-max(get(paste("vertex",i,sep=""))$Freq)
  vertex_gini[j]<-Gini(get(paste("vertex",i,sep=""))$Freq)
  cluster_max[j]<-max(table(get(paste("edges",i,sep=""))$igh_clone_id))
  clusters[j]<-sum(table(table(get(paste("edges",i,sep=""))$igh_clone_id)))
  num_reads_max_cluster[j]<-tail(table(table(get(paste("edges",i,sep=""))$igh_clone_id)),1)
  cluster_gini[j]<-Gini(table(get(paste("edges",i,sep=""))$igh_clone_id))
  j=j+1
}

id<-match(specimen, reads_clones_annot_AR_qc$specimen_id)
reads_clones_annot_AR_qc$cluster_gini[id]<-cluster_gini
reads_clones_annot_AR_qc$vertex_gini[id]<-vertex_gini
reads_clones_annot_AR_qc$vertex_max[id]<-vertex_max
reads_clones_annot_AR_qc$cluster_max[id]<-cluster_max
reads_clones_annot_AR_qc$num_reads_max_cluster[id]<-num_reads_max_cluster
reads_clones_annot_AR_qc$clusters[id]<-clusters
reads_clones_annot_AR_qc$clonal_expansion<-(reads_clones_annot_AR_qc$num_reads_max_cluster/reads_clones_annot_AR_qc$reads_gDNA)*100

reads_clones_annot_AR_qc$clin<-factor(reads_clones_annot_AR_qc$clin)
reads_clones_annot_AR_qc$clin<-relevel(reads_clones_annot_AR_qc$clin,ref="pre-AR")

tiff("plot_vertex_gini_cluster_gini_AR_gDNA.tiff",h=2000,w=1800,res=300)
qplot(cluster_gini, vertex_gini, data= reads_clones_annot_AR_qc, colour = clin) +
  scale_colour_manual(values=c("goldenrod","firebrick3"))
dev.off()

tiff("boxplot_vertex_gini_AR_gDNA.tiff",h=2000,w=1800,res=300)
ggplot(reads_clones_annot_AR_qc, aes(clin,vertex_gini,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Vertex Gini")
dev.off()
summary(glm(reads_clones_annot_AR_qc$vertex_gini ~ reads_clones_annot_AR_qc$clin))

tiff("boxplot_cluster_gini_AR_gDNA.tiff",h=2000,w=1800,res=300)
ggplot(reads_clones_annot_AR_qc, aes(clin,cluster_gini,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Cluster Gini")
dev.off()
summary(glm(reads_clones_annot_AR_qc$cluster_gini ~ reads_clones_annot_AR_qc$clin))

###Once is in network format I plot the network
##NP
specimen_preAR<- c("7_S80","7_S84","7_S86","7_S88","7_S90")
specimen_AR<- c("7_S81","7_S85","7_S87","7_S89","7_S91")

for(i in specimen_preAR) {
  print(i)
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("goldenrod")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}


for(i in specimen_AR) {
  print(i)
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("firebrick3")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}



############
### cDNA ###
############

data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

### longitudinal
data_qc_cDNA_long<-data_qc_cDNA[which(data_qc_cDNA$clin!="AR" & data_qc_cDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$clones_cDNA>1000),]

id<-match(data_qc_cDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_cDNA_long_qc<-data_qc_cDNA_long[which(is.na(id)==F),]
data_qc_cDNA_long_qc$specimen_label<-factor(data_qc_cDNA_long_qc$specimen_label)

data_qc_cDNA_long_qc$clone_CDR3<-paste(data_qc_cDNA_long_qc$igh_clone_id, data_qc_cDNA_long_qc$cdr3_seq_aa_q,sep="")
#write.table(data_qc_cDNA_long_qc,"/Users/Pinedasans/VDJ/Data/data_qc_cDNA_long_qc.txt",sep="\t",row.names = F)

###Save to apply the nucleotides-assembly-1.0.jar made by Mikel
specimen<-unique(data_qc_cDNA_long_qc$specimen_label)
for (i in specimen){
  print(i)
  data_qc_cDNA_long_qc_specimen<-data_qc_cDNA_long_qc[which(data_qc_cDNA_long_qc$specimen_label==i),]
  #write.table(data_qc_cDNA_long_qc_specimen,paste("/Users/Pinedasans/VDJ/Data/data_qc_cDNA_long_qc_",i,".txt",sep=""),sep="\t",row.names = F)
  df_specimen<-data_qc_cDNA_long_qc_specimen[,c("clone_CDR3","igh_clone_id")]
  groups <- group_by(df_specimen,clone_CDR3,igh_clone_id)
  edges<-unique(data.frame(groups))
  vertex<-data.frame(table(data_qc_cDNA_long_qc_specimen$clone_CDR3))
  write.table(edges,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/long_cDNA/edges",i,".txt",sep=""),sep="\t",row.names = F)
  write.table(vertex,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/long_cDNA/vertex",i,".txt",sep=""),sep="\t",row.names = F)
}


##Obtain the vertex and cluster Gini index
vertex_max<-NULL
vertex_gini<-NULL
cluster_max<-NULL
cluster_gini<-NULL
num_reads_max_cluster<-NULL
clusters<-NULL
j<-1
specimen<-unique(data_qc_cDNA_long_qc$specimen_label)
for (i in specimen){
  assign(paste("edges",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/long_cDNA/edges",i,".txt",sep=""),header=T))
  assign(paste("vertex",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/long_cDNA/vertex",i,".txt",sep=""),header=T))
  vertex_max[j]<-max(get(paste("vertex",i,sep=""))$Freq)
  vertex_gini[j]<-Gini(get(paste("vertex",i,sep=""))$Freq)
  cluster_max[j]<-max(table(get(paste("edges",i,sep=""))$igh_clone_id))
  clusters[j]<-sum(table(table(get(paste("edges",i,sep=""))$igh_clone_id)))
  num_reads_max_cluster[j]<-tail(table(table(get(paste("edges",i,sep=""))$igh_clone_id)),1)
  cluster_gini[j]<-Gini(table(get(paste("edges",i,sep=""))$igh_clone_id))
  j=j+1
}

 
id<-match(specimen, reads_clones_annot_Long_qc$specimen_id)
reads_clones_annot_Long_qc$cluster_gini[id]<-cluster_gini
reads_clones_annot_Long_qc$vertex_gini[id]<-vertex_gini
reads_clones_annot_Long_qc$vertex_max[id]<-vertex_max
reads_clones_annot_Long_qc$cluster_max[id]<-cluster_max
reads_clones_annot_Long_qc$num_reads_max_cluster[id]<-num_reads_max_cluster
reads_clones_annot_Long_qc$clusters[id]<-clusters
reads_clones_annot_Long_qc$clonal_expansion<-(reads_clones_annot_Long_qc$num_reads_max_cluster/reads_clones_annot_Long_qc$reads_cDNA)*100

reads_clones_annot_Long_qc$clin<-factor(reads_clones_annot_Long_qc$clin)

tiff("plot_vertex_gini_cluster_gini_gDNA.tiff",h=1600,w=2000,res=300)
reads_clones_annot_Long_qc$timeplot<-ifelse(reads_clones_annot_Long_qc$time2==0,1,
                                            ifelse(reads_clones_annot_Long_qc$time2==6,1.5,2))
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[reads_clones_annot_Long_qc$clin]
plot(cluster_gini, vertex_gini,cex=reads_clones_annot_Long_qc$timeplot,
     col = cols,pch=20,ylab = "Vextex Gini",xlab = "Cluster Gini")
legend("topleft",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
dev.off()


tiff("plot_vertex_gini_cluster_gini_cDNA.tiff",h=2000,w=1800,res=300)
reads_clones_annot_Long_qc_time0<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==0),]
p1<-qplot(cluster_gini, vertex_gini, data= reads_clones_annot_Long_qc_time0, colour = clin) +
  scale_colour_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) +
  scale_x_continuous(limit = c(0,0.25)) + scale_y_continuous(limit = c(0,0.75)) + labs(title="time 0")

reads_clones_annot_Long_qc_time6<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==6),]
p2<-qplot(cluster_gini, vertex_gini, data= reads_clones_annot_Long_qc_time6, colour = clin) +
  scale_colour_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) +
  scale_x_continuous(limit = c(0,0.25)) + scale_y_continuous(limit = c(0,0.75)) + labs(title="time 6")

reads_clones_annot_Long_qc_time24<-reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time2==24),]
p3<-qplot(cluster_gini, vertex_gini, data= reads_clones_annot_Long_qc_time24, colour = clin) +
  scale_colour_manual(values=c("chartreuse4","dodgerblue3","darkorange2")) +
  scale_x_continuous(limit = c(0,0.25)) + scale_y_continuous(limit = c(0,0.75)) + labs(title="time 24")
multiplot(p1,p2,p3)
dev.off()

anova(lm(reads_clones_annot_Long_qc$cluster_gini~reads_clones_annot_Long_qc$clin))
tiff("boxplot_vertexGini_clin_cDNA.tiff",h=2000,w=1800,res=300)
boxplot(reads_clones_annot_Long_qc$clusters~reads_clones_annot_Long_qc$clin,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="Degree of clonal expansion")
dev.off()


tiff("boxplot_vertex_gini_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time==0),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==0)]), 
                reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Vertex Gini")
p2 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time==6),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==6)]), 
                reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 6",x="Clin", y = "Vertex Gini")
p3 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time==24),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==24)]), 
                reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 24",x="Clin", y = "Vertex Gini")
multiplot(p1,p2,p3)
dev.off()

summary(glm(reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==0)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==0)]))
summary(glm(reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==6)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==6)]))
summary(glm(reads_clones_annot_Long_qc$vertex_gini[which(reads_clones_annot_Long_qc$time==24)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==24)]))


tiff("boxplot_cluster_gini_cDNA.tiff",h=2000,w=1800,res=300)
p1 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time==0),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==0)]), 
                reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==0)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 0",x="Clin", y = "Cluster Gini")
p2 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time==6),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==6)]), 
                reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==6)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 6",x="Clin", y = "Cluster Gini")
p3 = ggplot(reads_clones_annot_Long_qc[which(reads_clones_annot_Long_qc$time==24),], 
            aes(factor(reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==24)]), 
                reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==24)],fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + 
  labs(title="time 24",x="Clin", y = "Cluster Gini")
multiplot(p1,p2,p3)
dev.off()

summary(glm(reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==0)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==0)]))
summary(glm(reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==6)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==6)]))
summary(glm(reads_clones_annot_Long_qc$cluster_gini[which(reads_clones_annot_Long_qc$time==24)] ~ reads_clones_annot_Long_qc$clin[which(reads_clones_annot_Long_qc$time==24)]))


###Once is in network format I plot the network

##NP
files <- list.files("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/NP/")
for(i in files) {
  cat(i, "\n")
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/NP/",i, sep = ""))
  j<-substr(i,6,10)
  nodes<-get(paste("vertex",j,sep=""))
  net<-graph_from_data_frame(d=edges,vertices = nodes,directed=T)
  V(net)$size <- V(net)$Freq/200
  V(net)$color <- c("chartreuse4")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$arrow.width <- 0.1
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_",j,".tiff",sep=""),res=300,h=1000,w=1000)
  plot(net,vertex.label=NA,edge.curved=0.1,layout=layout_with_graphopt(net,niter = 800))
  dev.off()
}
#PNR
files <- list.files("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/PNR/")
for(i in files) {
  cat(i, "\n")
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/PNR/",i, sep = ""))
  j<-substr(i,6,10)
  nodes<-get(paste("vertex",j,sep=""))
  net<-graph_from_data_frame(d=edges,vertices = nodes,directed=T)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("dodgerblue3")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_",j,".tiff",sep=""),res=300,h=1500,w=1500)
  plot(net,vertex.label=NA,edge.curved=0.1,layout=layout_with_graphopt(net))
  dev.off()
}

##PR
files <- list.files("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/PR/")
for(i in files) {
  cat(i, "\n")
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/PR/",i, sep = ""))
  j<-substr(i,6,10)
  nodes<-get(paste("vertex",j,sep=""))
  net<-graph_from_data_frame(d=edges,vertices = nodes,directed=T)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("darkorange2")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_",j,".tiff",sep=""),res=300,h=1500,w=1500)
  plot(net,vertex.label=NA,edge.curved=0.1,layout=layout_with_graphopt(net))
  dev.off()
}

##############
### AR #######
#############
data_qc_cDNA_AR<-data_qc_cDNA[which(data_qc_cDNA$clin=="AR" | data_qc_cDNA$clin=="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
reads_clones_annot_AR_qc<-reads_clones_annot_AR[which(reads_clones_annot_AR$reads_cDNA>100),]

id<-match(data_qc_cDNA_AR$specimen_label,reads_clones_annot_AR_qc$specimen_id)
data_qc_cDNA_AR_qc<-data_qc_cDNA_AR[which(is.na(id)==F),]
data_qc_cDNA_AR_qc$specimen_label<-factor(data_qc_cDNA_AR_qc$specimen_label)

data_qc_cDNA_AR_qc$clone_CDR3<-paste(data_qc_cDNA_AR_qc$igh_clone_id, data_qc_cDNA_AR_qc$cdr3_seq_aa_q,sep="")


###Save to apply the nucleotides-assembly-1.0.jar made by Mikel
specimen<-unique(data_qc_cDNA_AR_qc$specimen_label)
for (i in specimen){
  print(i)
  data_qc_cDNA_AR_qc_specimen<-data_qc_cDNA_AR_qc[which(data_qc_cDNA_AR_qc$specimen_label==i),]
  df_specimen<-data_qc_cDNA_AR_qc_specimen[,c("clone_CDR3","igh_clone_id")]
  groups <- group_by(df_specimen,clone_CDR3,igh_clone_id)
  edges<-unique(data.frame(groups))
  vertex<-data.frame(table(data_qc_cDNA_AR_qc_specimen$clone_CDR3))
  write.table(edges,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/AR_cDNA/edges",i,".txt",sep=""),sep="\t",row.names = F)
  write.table(vertex,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_cDNA//vertex",i,".txt",sep=""),sep="\t",row.names = F)
}

##Obtain the vertex and cluster Gini index
vertex_max<-NULL
vertex_gini<-NULL
cluster_max<-NULL
cluster_gini<-NULL
num_reads_max_cluster<-NULL
clusters<-NULL
j<-1
specimen<-unique(data_qc_cDNA_AR_qc$specimen_label)
for (i in specimen){
  assign(paste("edges",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/AR_cDNA/edges",i,".txt",sep=""),header=T))
  assign(paste("vertex",i,sep=""),read.table(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_AR_cDNA/vertex",i,".txt",sep=""),header=T))
  vertex_max[j]<-max(get(paste("vertex",i,sep=""))$Freq)
  vertex_gini[j]<-Gini(get(paste("vertex",i,sep=""))$Freq)
  cluster_max[j]<-max(table(get(paste("edges",i,sep=""))$igh_clone_id))
  clusters[j]<-sum(table(table(get(paste("edges",i,sep=""))$igh_clone_id)))
  num_reads_max_cluster[j]<-tail(table(table(get(paste("edges",i,sep=""))$igh_clone_id)),1)
  cluster_gini[j]<-Gini(table(get(paste("edges",i,sep=""))$igh_clone_id))
  j=j+1
}

id<-match(specimen, reads_clones_annot_AR_qc$specimen_id)
reads_clones_annot_AR_qc$cluster_gini[id]<-cluster_gini
reads_clones_annot_AR_qc$vertex_gini[id]<-vertex_gini
reads_clones_annot_AR_qc$vertex_max[id]<-vertex_max
reads_clones_annot_AR_qc$cluster_max[id]<-cluster_max
reads_clones_annot_AR_qc$num_reads_max_cluster[id]<-num_reads_max_cluster
reads_clones_annot_AR_qc$clusters[id]<-clusters
reads_clones_annot_AR_qc$clonal_expansion<-(reads_clones_annot_AR_qc$num_reads_max_cluster/reads_clones_annot_AR_qc$reads_cDNA)*100

reads_clones_annot_AR_qc$clin<-factor(reads_clones_annot_AR_qc$clin)
reads_clones_annot_AR_qc$clin<-relevel(reads_clones_annot_AR_qc$clin,ref="pre-AR")

tiff("plot_vertex_gini_cluster_gini_AR_cDNA.tiff",h=2000,w=1800,res=300)
qplot(cluster_gini, vertex_gini, data= reads_clones_annot_AR_qc, colour = clin) +
  scale_colour_manual(values=c("goldenrod","firebrick3"))
dev.off()

tiff("boxplot_vertex_gini_AR_cDNA.tiff",h=2000,w=1800,res=300)
ggplot(reads_clones_annot_AR_qc, aes(clin,vertex_gini,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Vertex Gini")
dev.off()
summary(glm(reads_clones_annot_AR_qc$vertex_gini ~ reads_clones_annot_AR_qc$clin))

tiff("boxplot_cluster_gini_AR_cDNA.tiff",h=2000,w=1800,res=300)
ggplot(reads_clones_annot_AR_qc, aes(clin,cluster_gini,fill=clin)) + 
  geom_boxplot() + scale_fill_manual(values=c("goldenrod","firebrick3")) + labs(x="Clin", y = "Cluster Gini")
dev.off()
summary(glm(reads_clones_annot_AR_qc$cluster_gini ~ reads_clones_annot_AR_qc$clin))















clone_df<-data.frame(table(data_qc_cDNA_long_qc$igh_clone_id,data_qc_cDNA_long_qc$specimen_label))
colnames(clone_df)<-c("clone","specimen","count")
clone_df$clin = reads_clones_annot_Long_qc[clone_df$specimen,1]
clone_df$time = reads_clones_annot_Long_qc[clone_df$specimen,4]
clone_df$sample = reads_clones_annot_Long_qc[clone_df$specimen,2]

clone_df_noceros = clone_df[which(clone_df$count!=0),]


g1<-ggplot(clone_df_noceros[which(clone_df_noceros$clin=="NP"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clone Expansion")

g2<-ggplot(clone_df_noceros[which(clone_df_noceros$clin=="PNR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clone Expansion")

g3<-ggplot(clone_df_noceros[which(clone_df_noceros$clin=="PR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clone Expansion")

tiff("clone_expansion_cDNA_long.tiff",res=300,h=1700,w=2700)
multiplot(g1,g2,g3)
dev.off()




v_gene_df<-data.frame(table(data_qc_cDNA_long_qc$v_gene,data_qc_cDNA_long_qc$specimen_label))
colnames(v_gene_df)<-c("v_gene","specimen","count")
v_gene_df$clin = reads_clones_annot_Long_qc[v_gene_df$specimen,1] 
v_gene_df$time = reads_clones_annot_Long_qc[v_gene_df$specimen,4] 
v_gene_df$sample = reads_clones_annot_Long_qc[v_gene_df$specimen,2]
v_gene_df$clones = reads_clones_annot_Long_qc[v_gene_df$specimen,42]


g1<-ggplot(v_gene_df[which(v_gene_df$clin=="NP"),], aes(x=time,y=count/clones,group=v_gene,fill=v_gene)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  scale_y_continuous(limit = c(0,300000)) + facet_grid(clin ~ sample) + labs(x = "time", y = "V gene usage")

g2<-ggplot(v_gene_df[which(v_gene_df$clin=="PNR"),], aes(x=time,y=count/clones,group=v_gene,fill=v_gene)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  scale_y_continuous(limit = c(0,300000)) + facet_grid(clin ~ sample) + labs(x = "time", y = "V gene usage")

g3<-ggplot(v_gene_df[which(v_gene_df$clin=="PR"),], aes(x=time,y=count/clones,group=v_gene,fill=v_gene)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  scale_y_continuous(limit = c(0,300000)) + facet_grid(clin ~ sample) + labs(x = "time", y = "V gene usage")

tiff("v_gene_usage_cDNA_long.tiff",res=300,h=1700,w=2700)
multiplot(g1,g2,g3)
dev.off()


### AR
data_qc_cDNA_AR<-data_qc_cDNA[which(data_qc_cDNA$clin=="AR" | data_qc_cDNA$clin=="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
reads_clones_annot_AR_qc<-reads_clones_annot_AR[which(reads_clones_annot_AR$reads_cDNA>100),]

id<-match(data_qc_cDNA_AR$specimen_label,reads_clones_annot_AR_qc$specimen_id)
data_qc_cDNA_AR_qc<-data_qc_cDNA_AR[which(is.na(id)==F),]
data_qc_cDNA_AR_qc$specimen_label<-factor(data_qc_cDNA_AR_qc$specimen_label)

v_gene_df<-data.frame(table(data_qc_cDNA_AR_qc$v_gene,data_qc_cDNA_AR_qc$specimen_label))
colnames(v_gene_df)<-c("v_gene","specimen","count")
v_gene_df$clin = reads_clones_annot_AR_qc[v_gene_df$specimen,1] 
v_gene_df$time = reads_clones_annot_AR_qc[v_gene_df$specimen,4] 
v_gene_df$sample = reads_clones_annot_AR_qc[v_gene_df$specimen,2]
v_gene_df$clin<-relevel(v_gene_df$clin,ref="pre-AR")

tiff("v_gene_usage_cDNA_AR.tiff",res=300,h=1700,w=2700)
ggplot(v_gene_df, aes(x=clin,y=count,group=v_gene,fill=v_gene)) +
  geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  facet_grid(~ sample) + labs(x = "time", y = "V gene usage")
dev.off()


###### V-gene usage analysis ####
data_qc_gDNA$specimen_label<-factor(data_qc_gDNA$specimen_label)
v_gene_matrix<-table(data_qc_gDNA$v_gene,data_qc_gDNA$specimen_label)
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
id<-match(reads_clones_annot_Long$specimen_id,colnames(v_gene_matrix))
v_gene_matrix_long<-v_gene_matrix[,na.omit(id)]
id2<-match(colnames(v_gene_matrix_long),reads_clones_annot_Long$specimen_id)

p.value<-NULL
for (i in 1:nrow(v_gene_matrix_long)){
  result<-anova(lm(v_gene_matrix_long[i,]~reads_clones_annot_Long$clin[id2]))
  p.value[i]<-result$`Pr(>F)`[1]
}
p.adj<-p.adjust(p.value,method = "fdr")
v_gene_matrix_long_sign<-v_gene_matrix_long[which(p.adj<0.1),]

annot_col<-data.frame(reads_clones_annot_Long$clin[id2])
rownames(annot_col)<-colnames(v_gene_matrix_long_sign)
colnames(annot_col)<-"clin"
pheatmap(v_gene_matrix_long_sign,scale = "row",border_color = NA,annotation_col = annot_col )




