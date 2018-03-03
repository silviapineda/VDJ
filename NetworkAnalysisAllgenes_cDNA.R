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

setwd("/Users/Pinedasans/VDJ/ResultsAllClones/")
load("~/VDJ/Data/VDJ_clonesAllmerged.Rdata")


##########################
### cDNA isotype IGH- ###
#########################
data_cDNA<-data_merge[which(data_merge$amplification_template=="cDNA"),]
data_cDNA_IGHA<-data_merge[which(data_merge$isotype=="IGHA"),]
data_cDNA_long<-data_cDNA_IGHA[which(data_cDNA_IGHA$clin!="AR" & data_cDNA_IGHA$clin!="pre-AR"),]
reads_clones_annot_long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_long_cDNA<-reads_clones_annot_long[which(reads_clones_annot_long$clones_IGHA>100),]
id<-match(data_cDNA_long$specimen_label,reads_clones_annot_long_cDNA$specimen_id)
data_cDNA_long_qc<-data_cDNA_long[which(is.na(id)==F),]
data_cDNA_long_qc$specimen_label<-factor(data_cDNA_long_qc$specimen_label)

#################################
###### Network Analysis  ########
#################################
##First: obtain the vertex and edges
data_cDNA_long_qc$CloneId_CDR3<-paste0(data_cDNA_long_qc[,c("V_J_lenghCDR3_CloneId")],data_cDNA_long_qc[,c("cdr3_seq_aa_q")])
specimen<-unique(data_cDNA_long_qc$specimen_label)
for (i in specimen){
  print(i)
  data_cDNA_long_qc_specimen<-data_cDNA_long_qc[which(data_cDNA_long_qc$specimen_label==i),]
  df_specimen<-data_cDNA_long_qc_specimen[,c("CloneId_CDR3","V_J_lenghCDR3_CloneId")]
  groups <- group_by(df_specimen,CloneId_CDR3,V_J_lenghCDR3_CloneId)
  assign(paste0("edges",i),unique(data.frame(groups)))
  df_vertex<-data.frame(table(data_cDNA_long_qc_specimen$CloneId_CDR3))
  assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
  #write.table(get(paste0("edges",i)),paste0("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHA/edges",i,".txt"),sep="\t",row.names = F)
  #write.table(get(paste0("vertex",i)),paste0("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHA/vertex",i,".txt"),sep="\t",row.names = F)
}

##2. Apply the nucleotides-assembly-1.0.jar made by Mikel

##3. Obtain the vertex and cluster Gini index
vertex_max<-NULL
vertex_gini<-NULL
cluster_max<-NULL
cluster_gini<-NULL
num_reads_max_cluster<-NULL
clusters<-NULL
j<-1
specimen<-unique(data_cDNA_long_qc$specimen_label)
for (i in specimen){
  vertex_max[j]<-max(get(paste("vertex",i,sep=""))$Freq)
  vertex_gini[j]<-Gini(get(paste("vertex",i,sep=""))$Freq)
  cluster_max[j]<-max(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId))
  clusters[j]<-sum(table(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId)))
  num_reads_max_cluster[j]<-tail(table(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId)),1)
  cluster_gini[j]<-Gini(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId))
  j=j+1
}
names(cluster_gini)<-specimen

id<-match(specimen, reads_clones_annot_long_cDNA$specimen_id)
reads_clones_annot_long_cDNA$cluster_gini[id]<-cluster_gini
reads_clones_annot_long_cDNA$vertex_gini[id]<-vertex_gini
reads_clones_annot_long_cDNA$vertex_max[id]<-vertex_max
reads_clones_annot_long_cDNA$cluster_max[id]<-cluster_max
reads_clones_annot_long_cDNA$num_reads_max_cluster[id]<-num_reads_max_cluster
reads_clones_annot_long_cDNA$clusters[id]<-clusters
reads_clones_annot_long_cDNA$clonal_expansion<-(reads_clones_annot_long_cDNA$num_reads_max_cluster/reads_clones_annot_long_cDNA$reads_cDNA)*100

reads_clones_annot_long_cDNA$clin<-factor(reads_clones_annot_long_cDNA$clin)

##change the 12 for the 6                                                                             
##Delete sample M154-S001
reads_clones_annot_long_cDNA<-reads_clones_annot_long_cDNA[which(reads_clones_annot_long_cDNA$specimen_id!="M154-S001"),]
reads_clones_annot_long_cDNA$time<-replace(reads_clones_annot_long_cDNA$time,reads_clones_annot_long_cDNA$time==12,6)
##Change the 13 for 24
reads_clones_annot_long_cDNA$time<-replace(reads_clones_annot_long_cDNA$time,reads_clones_annot_long_cDNA$time==13,24)

#Time 9 delete and time 32
reads_clones_annot_long_cDNA<-reads_clones_annot_long_cDNA[which(reads_clones_annot_long_cDNA$time!=9),]
reads_clones_annot_long_cDNA<-reads_clones_annot_long_cDNA[which(reads_clones_annot_long_cDNA$time!=32),]

reads_clones_annot_long_cDNA$timeplot<-ifelse(reads_clones_annot_long_cDNA$time==0,1,
                               ifelse(reads_clones_annot_long_cDNA$time==6,1.5,2))
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[reads_clones_annot_long_cDNA$clin]

tiff("network_vertex_cluster_gini_cDNA_IGHA.tiff",h=2000,w=2000,res=300)
par(fig=c(0,0.8,0,0.8))
plot(reads_clones_annot_long_cDNA$cluster_gini, reads_clones_annot_long_cDNA$vertex_gini,
     cex=reads_clones_annot_long_cDNA$timeplot,
        col = cols,pch=20,ylab = "Vextex Gini",xlab = "Cluster Gini")
legend("bottomright",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)

reads_clones_annot_Long_qc_time24<-reads_clones_annot_long_cDNA[which(reads_clones_annot_long_cDNA$time==24),]
par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(lm(reads_clones_annot_Long_qc_time24$cluster_gini~reads_clones_annot_Long_qc_time24$clin))
boxplot(reads_clones_annot_Long_qc_time24$cluster_gini~reads_clones_annot_Long_qc_time24$clin,
        col=c("chartreuse4", "dodgerblue3","darkorange2"), horizontal=TRUE, axes=FALSE)

par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(lm(reads_clones_annot_Long_qc_time24$vertex_gini~reads_clones_annot_Long_qc_time24$clin))
boxplot(reads_clones_annot_Long_qc_time24$vertex_gini~reads_clones_annot_Long_qc_time24$clin,
        col=c("chartreuse4", "dodgerblue3","darkorange2"),axes=FALSE)
dev.off()

##Vertex gini 
summary(glm(reads_clones_annot_long_cDNA$vertex_gini[which(reads_clones_annot_long_cDNA$time==0)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==0)]))
summary(glm(reads_clones_annot_long_cDNA$vertex_gini[which(reads_clones_annot_long_cDNA$time==6)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==6)]))
summary(glm(reads_clones_annot_long_cDNA$vertex_gini[which(reads_clones_annot_long_cDNA$time==24)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==24)]))

summary(glm(reads_clones_annot_long_cDNA$cluster_gini[which(reads_clones_annot_long_cDNA$time==0)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==0)]))
summary(glm(reads_clones_annot_long_cDNA$cluster_gini[which(reads_clones_annot_long_cDNA$time==6)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==6)]))
summary(glm(reads_clones_annot_long_cDNA$cluster_gini[which(reads_clones_annot_long_cDNA$time==24)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==24)]))


##############################
### cDNA naive and memory ###
#############################
data_cDNA<-data_merge[which(data_merge$amplification_template=="cDNA"),]
data_cDNA_memory<-data_merge[which(data_merge$IGHM_naive_memory=="memory"),]
data_cDNA_long<-data_cDNA_memory[which(data_cDNA_memory$clin!="AR" & data_cDNA_memory$clin!="pre-AR"),]
reads_clones_annot_long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_long_cDNA<-reads_clones_annot_long[which(reads_clones_annot_long$clones_memory>100),]
id<-match(data_cDNA_long$specimen_label,reads_clones_annot_long_cDNA$specimen_id)
data_cDNA_long_qc<-data_cDNA_long[which(is.na(id)==F),]
data_cDNA_long_qc$specimen_label<-factor(data_cDNA_long_qc$specimen_label)

#################################
###### Network Analysis  ########
#################################
##First: obtain the vertex and edges
data_cDNA_long_qc$CloneId_CDR3<-paste0(data_cDNA_long_qc[,c("V_J_lenghCDR3_CloneId")],data_cDNA_long_qc[,c("cdr3_seq_aa_q")])
specimen<-unique(data_cDNA_long_qc$specimen_label)
for (i in specimen){
  print(i)
  data_cDNA_long_qc_specimen<-data_cDNA_long_qc[which(data_cDNA_long_qc$specimen_label==i),]
  df_specimen<-data_cDNA_long_qc_specimen[,c("CloneId_CDR3","V_J_lenghCDR3_CloneId")]
  groups <- group_by(df_specimen,CloneId_CDR3,V_J_lenghCDR3_CloneId)
  assign(paste0("edges",i),unique(data.frame(groups)))
  df_vertex<-data.frame(table(data_cDNA_long_qc_specimen$CloneId_CDR3))
  assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
  write.table(get(paste0("edges",i)),paste0("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_memory/edges",i,".txt"),sep="\t",row.names = F)
  write.table(get(paste0("vertex",i)),paste0("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_memory/vertex",i,".txt"),sep="\t",row.names = F)
}

##2. Apply the nucleotides-assembly-1.0.jar made by Mikel

##3. Obtain the vertex and cluster Gini index
vertex_max<-NULL
vertex_gini<-NULL
cluster_max<-NULL
cluster_gini<-NULL
num_reads_max_cluster<-NULL
clusters<-NULL
j<-1
specimen<-unique(data_cDNA_long_qc$specimen_label)
for (i in specimen){
  vertex_max[j]<-max(get(paste("vertex",i,sep=""))$Freq)
  vertex_gini[j]<-Gini(get(paste("vertex",i,sep=""))$Freq)
  cluster_max[j]<-max(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId))
  clusters[j]<-sum(table(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId)))
  num_reads_max_cluster[j]<-tail(table(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId)),1)
  cluster_gini[j]<-Gini(table(get(paste("edges",i,sep=""))$V_J_lenghCDR3_CloneId))
  j=j+1
}
names(cluster_gini)<-specimen

id<-match(specimen, reads_clones_annot_long_cDNA$specimen_id)
reads_clones_annot_long_cDNA$cluster_gini[id]<-cluster_gini
reads_clones_annot_long_cDNA$vertex_gini[id]<-vertex_gini
reads_clones_annot_long_cDNA$vertex_max[id]<-vertex_max
reads_clones_annot_long_cDNA$cluster_max[id]<-cluster_max
reads_clones_annot_long_cDNA$num_reads_max_cluster[id]<-num_reads_max_cluster
reads_clones_annot_long_cDNA$clusters[id]<-clusters
reads_clones_annot_long_cDNA$clonal_expansion<-(reads_clones_annot_long_cDNA$num_reads_max_cluster/reads_clones_annot_long_cDNA$reads_cDNA)*100

reads_clones_annot_long_cDNA$clin<-factor(reads_clones_annot_long_cDNA$clin)

##change the 12 for the 6                                                                             
reads_clones_annot_long_cDNA$time<-replace(reads_clones_annot_long_cDNA$time,reads_clones_annot_long_cDNA$time==12,6)
##Change the 13 for 24
reads_clones_annot_long_cDNA$time<-replace(reads_clones_annot_long_cDNA$time,reads_clones_annot_long_cDNA$time==13,24)

#Time 9 delete and time 32
reads_clones_annot_long_cDNA<-reads_clones_annot_long_cDNA[which(reads_clones_annot_long_cDNA$time!=9),]
reads_clones_annot_long_cDNA<-reads_clones_annot_long_cDNA[which(reads_clones_annot_long_cDNA$time!=32),]

reads_clones_annot_long_cDNA$timeplot<-ifelse(reads_clones_annot_long_cDNA$time==0,1,
                                              ifelse(reads_clones_annot_long_cDNA$time==6,1.5,2))
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[reads_clones_annot_long_cDNA$clin]

tiff("network_vertex_cluster_gini_cDNA_memory.tiff",h=2000,w=2000,res=300)
par(fig=c(0,0.8,0,0.8))
plot(reads_clones_annot_long_cDNA$cluster_gini, reads_clones_annot_long_cDNA$vertex_gini,
     cex=reads_clones_annot_long_cDNA$timeplot,
     col = cols,pch=20,ylab = "Vextex Gini",xlab = "Cluster Gini")
#legend("bottomright",legend=c("NP","PNR","PR","t0","t6","t24"), 
#      col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
#     pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)

reads_clones_annot_Long_qc_time24<-reads_clones_annot_long_cDNA[which(reads_clones_annot_long_cDNA$time==24),]
par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(lm(reads_clones_annot_Long_qc_time24$cluster_gini~reads_clones_annot_Long_qc_time24$clin))
boxplot(reads_clones_annot_Long_qc_time24$cluster_gini~reads_clones_annot_Long_qc_time24$clin,
        col=c("chartreuse4", "dodgerblue3","darkorange2"), horizontal=TRUE, axes=FALSE)

par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(lm(reads_clones_annot_Long_qc_time24$vertex_gini~reads_clones_annot_Long_qc_time24$clin))
boxplot(reads_clones_annot_Long_qc_time24$vertex_gini~reads_clones_annot_Long_qc_time24$clin,
        col=c("chartreuse4", "dodgerblue3","darkorange2"),axes=FALSE)
dev.off()

##Vertex gini 
summary(glm(reads_clones_annot_long_cDNA$vertex_gini[which(reads_clones_annot_long_cDNA$time==0)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==0)]))
summary(glm(reads_clones_annot_long_cDNA$vertex_gini[which(reads_clones_annot_long_cDNA$time==6)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==6)]))
summary(glm(reads_clones_annot_long_cDNA$vertex_gini[which(reads_clones_annot_long_cDNA$time==24)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==24)]))

summary(glm(reads_clones_annot_long_cDNA$cluster_gini[which(reads_clones_annot_long_cDNA$time==0)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==0)]))
summary(glm(reads_clones_annot_long_cDNA$cluster_gini[which(reads_clones_annot_long_cDNA$time==6)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==6)]))
summary(glm(reads_clones_annot_long_cDNA$cluster_gini[which(reads_clones_annot_long_cDNA$time==24)] ~ reads_clones_annot_long_cDNA$clin[which(reads_clones_annot_long_cDNA$time==24)]))



###Once is in network format I plot the network for IGHM which is the only one significant
##NP
specimen_NP<-as.character(reads_clones_annot_long_cDNA$specimen_id[which(reads_clones_annot_long_cDNA$clin=="NP")])
specimen_PNR<-reads_clones_annot_long_cDNA$specimen_id[which(reads_clones_annot_long_cDNA$clin=="PNR")]
specimen_PR<-reads_clones_annot_long_cDNA$specimen_id[which(reads_clones_annot_long_cDNA$clin=="PR")]

for(i in specimen_NP) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("chartreuse4")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/network_cDNA_IGHM",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}

for(i in specimen_PNR) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("dodgerblue3")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/network_cDNA_IGHM",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}

for(i in specimen_PR) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("darkorange2")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_cDNA_IGHM/network_cDNA_IGHM",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}
