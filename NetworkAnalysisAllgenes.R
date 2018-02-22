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
library(gridExtra)
library(grid)
library(lattice)

setwd("/Users/Pinedasans/VDJ/ResultsAllClones/")
load("~/VDJ/Data/VDJ_clonesAllmerged.Rdata")


############
### gDNA ###
############
data_gDNA<-data_merge[which(data_merge$amplification_template=="gDNA"),]
data_gDNA_long<-data_gDNA[which(data_gDNA$clin!="AR" & data_gDNA$clin!="pre-AR"),]
reads_clones_annot_long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_long_gDNA<-reads_clones_annot_long[which(reads_clones_annot_long$clones_gDNA>=100),]
id<-match(data_gDNA_long$specimen_label,reads_clones_annot_long_gDNA$specimen_id)
data_gDNA_long_qc<-data_gDNA_long[which(is.na(id)==F),]
data_gDNA_long_qc$specimen_label<-factor(data_gDNA_long_qc$specimen_label)

#################################
###### Network Analysis  ########
#################################
##First: obtain the vertex and edges
data_gDNA_long_qc$CloneId_CDR3<-paste0(data_gDNA_long_qc[,c("V_J_lenghCDR3_CloneId")],data_gDNA_long_qc[,c("cdr3_seq_aa_q")])
specimen<-unique(data_gDNA_long_qc$specimen_label)
for (i in specimen){
  print(i)
  data_gDNA_long_qc_specimen<-data_gDNA_long_qc[which(data_gDNA_long_qc$specimen_label==i),]
  df_specimen<-data_gDNA_long_qc_specimen[,c("CloneId_CDR3","V_J_lenghCDR3_CloneId")]
  groups <- group_by(df_specimen,CloneId_CDR3,V_J_lenghCDR3_CloneId)
  assign(paste0("edges",i),unique(data.frame(groups)))
  df_vertex<-data.frame(table(data_gDNA_long_qc_specimen$CloneId_CDR3))
  assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
  #write.table(get(paste0("edges",i)),paste0("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/edges",i,".txt"),sep="\t",row.names = F)
  #write.table(get(paste0("vertex",i)),paste0("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/vertex",i,".txt"),sep="\t",row.names = F)
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
specimen<-unique(data_gDNA_long_qc$specimen_label)
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

id<-match(specimen, reads_clones_annot_long_gDNA$specimen_id)
reads_clones_annot_long_gDNA$cluster_gini[id]<-cluster_gini
reads_clones_annot_long_gDNA$vertex_gini[id]<-vertex_gini
reads_clones_annot_long_gDNA$vertex_max[id]<-vertex_max
reads_clones_annot_long_gDNA$cluster_max[id]<-cluster_max
reads_clones_annot_long_gDNA$num_reads_max_cluster[id]<-num_reads_max_cluster
reads_clones_annot_long_gDNA$clusters[id]<-clusters
reads_clones_annot_long_gDNA$clonal_expansion<-(reads_clones_annot_long_gDNA$num_reads_max_cluster/reads_clones_annot_long_gDNA$reads_gDNA)*100


reads_clones_annot_long_gDNA$clin<-factor(reads_clones_annot_long_gDNA$clin)

###Delete sample 8 from the analysis
reads_clones_annot_long_gDNA<-reads_clones_annot_long_gDNA[which(reads_clones_annot_long_gDNA$Individual.id!="Individual8"),]


##change the 12 for the 6                                                                             
reads_clones_annot_long_gDNA$time<-replace(reads_clones_annot_long_gDNA$time,reads_clones_annot_long_gDNA$time==12,6)
##Change the 13 for 24
reads_clones_annot_long_gDNA$time<-replace(reads_clones_annot_long_gDNA$time,reads_clones_annot_long_gDNA$time==13,24)

#Time 9 delete and time 32
reads_clones_annot_long_gDNA<-reads_clones_annot_long_gDNA[which(reads_clones_annot_long_gDNA$time!=9),]
reads_clones_annot_long_gDNA<-reads_clones_annot_long_gDNA[which(reads_clones_annot_long_gDNA$time!=32),]

reads_clones_annot_long_gDNA$timeplot<-ifelse(reads_clones_annot_long_gDNA$time==0,1,
                               ifelse(reads_clones_annot_long_gDNA$time==6,1.5,2))
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
cols = COLOR[reads_clones_annot_long_gDNA$clin]

tiff("network_vertex_cluster_gini.tiff",h=2000,w=2000,res=300)
par(fig=c(0,0.8,0,0.8))
plot(reads_clones_annot_long_gDNA$cluster_gini, reads_clones_annot_long_gDNA$vertex_gini,
     cex=reads_clones_annot_long_gDNA$timeplot,
        col = cols,pch=20,ylab = "Gini (Vextex)",xlab = "Gini (Cluster)",xlim=c(0.0,0.09),ylim=c(0.2,0.55))
legend("bottomright",legend=c("NP","PNR","PR","t0","t6","t24"), 
       col=c("chartreuse4", "dodgerblue3","darkorange2","black","black","black"), 
       pch=20,cex=c(1.2),pt.cex=c(1.6,1.6,1.6,1,1.5,2),ncol=2)
reads_clones_annot_Long_qc_time24<-reads_clones_annot_long_gDNA[which(reads_clones_annot_long_gDNA$time==24),]

par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(lm(reads_clones_annot_Long_qc_time24$cluster_gini~reads_clones_annot_Long_qc_time24$clin))
boxplot(reads_clones_annot_Long_qc_time24$cluster_gini~reads_clones_annot_Long_qc_time24$clin,
        col=c("chartreuse4", "dodgerblue3","darkorange2"), horizontal=TRUE, axes=FALSE,
        ylim=c(0.0,0.09))

par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(lm(reads_clones_annot_Long_qc_time24$vertex_gini~reads_clones_annot_Long_qc_time24$clin))
boxplot(reads_clones_annot_Long_qc_time24$vertex_gini~reads_clones_annot_Long_qc_time24$clin,
        col=c("chartreuse4", "dodgerblue3","darkorange2"),axes=FALSE,ylim=c(0.2,0.55))
dev.off()

scatter <- ggplot(data=reads_clones_annot_long_gDNA,aes(x=cluster_gini, y=vertex_gini)) + geom_point(aes(color=clin)) + 
  scale_color_manual(values =c("chartreuse4", "dodgerblue3","darkorange2")) +
  theme(legend.position = "none")

plot_top <- ggplot(data=reads_clones_annot_Long_qc_time24, aes(x=cluster_gini, fill=clin)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c("chartreuse4", "dodgerblue3","darkorange2")) + 
  theme(legend.position = "none")

plot_right <- ggplot(data=reads_clones_annot_Long_qc_time24, aes(x=vertex_gini, fill=clin)) + 
  geom_density(alpha=.5) + coord_flip() +
  scale_fill_manual(values = c("chartreuse4", "dodgerblue3","darkorange2")) + 
  theme(legend.position = "none")
  
empty <- ggplot()+geom_point(aes(1,1), color="white") +
    theme(                              
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  
# Arrange the plots together
tiff("network_vertex_cluster_gini_type2.tiff",h=2000,w=2500,res=300)
grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

##Vertex gini 
summary(glm(reads_clones_annot_long_gDNA$vertex_gini[which(reads_clones_annot_long_gDNA$time==0)] ~ reads_clones_annot_long_gDNA$clin[which(reads_clones_annot_long_gDNA$time==0)]))
summary(glm(reads_clones_annot_long_gDNA$vertex_gini[which(reads_clones_annot_long_gDNA$time==6)] ~ reads_clones_annot_long_gDNA$clin[which(reads_clones_annot_long_gDNA$time==6)]))
summary(glm(reads_clones_annot_long_gDNA$vertex_gini[which(reads_clones_annot_long_gDNA$time==24)] ~ reads_clones_annot_long_gDNA$clin[which(reads_clones_annot_long_gDNA$time==24)]))

summary(glm(reads_clones_annot_long_gDNA$cluster_gini[which(reads_clones_annot_long_gDNA$time==0)] ~ reads_clones_annot_long_gDNA$clin[which(reads_clones_annot_long_gDNA$time==0)]))
summary(glm(reads_clones_annot_long_gDNA$cluster_gini[which(reads_clones_annot_long_gDNA$time==6)] ~ reads_clones_annot_long_gDNA$clin[which(reads_clones_annot_long_gDNA$time==6)]))
summary(glm(reads_clones_annot_long_gDNA$cluster_gini[which(reads_clones_annot_long_gDNA$time==24)] ~ reads_clones_annot_long_gDNA$clin[which(reads_clones_annot_long_gDNA$time==24)]))

##Cluster gini
fm_null <- lmer(reads_clones_annot_long_gDNA$cluster_gini ~ clin + time + (time | Sample_id),data=reads_clones_annot_long_gDNA,REML = F)
fm_full <- lmer(reads_clones_annot_long_gDNA$cluster_gini ~  clin*time + (time | Sample_id) ,data=reads_clones_annot_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_clusterGini_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = reads_clones_annot_long_gDNA$cluster_gini, colour = clin)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=3) +
  geom_smooth(method="lm",size=1.5) +
  labs(x = "time (months)",y = "Clones") 

print(p)
dev.off()

fm_null <- lmer(reads_clones_annot_long_gDNA$vertex_gini ~ clin + time + (time | Sample_id),data=reads_clones_annot_long_gDNA,REML = F)
fm_full <- lmer(reads_clones_annot_long_gDNA$vertex_gini ~  clin*time + (time | Sample_id) ,data=reads_clones_annot_long_gDNA,REML = F)
anova(fm_full, fm_null) 

tiff("plot_lmer_vertexGini_gDNA.tiff",h=1700,w=2000,res=300)
p <- ggplot(fm_full, aes(x = time, y = reads_clones_annot_long_gDNA$vertex_gini, colour = clin)) +
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
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("chartreuse4")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}

for(i in specimen_PNR) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("dodgerblue3")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}

for(i in specimen_PR) {
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/edges",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/long_gDNA/vertex",i, ".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("darkorange2")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("/Users/Pinedasans/VDJ/ResultsAllClones/network_data/network_gDNA",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}



