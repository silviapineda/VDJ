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
### Date: June, 2017
############################################################################################
library(ggplot2)
library(untb)
library(lme4)
library("caroline")
library("effects")
library("pheatmap")
library("igraph")

setwd("/Users/Pinedasans/VDJ/Results/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")



############
### gDNA ###
############
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
#data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

### longitudinal
data_qc_gDNA_long<-data_qc_gDNA[which(data_qc_gDNA$clin!="AR" & data_qc_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_gDNA>100),]

id<-match(data_qc_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_gDNA_long_qc<-data_qc_gDNA_long[which(is.na(id)==F),]
data_qc_gDNA_long_qc$specimen_label<-factor(data_qc_gDNA_long_qc$specimen_label)

data_qc_gDNA_long_qc$clone_CDR3<-paste(data_qc_gDNA_long_qc$igh_clone_id, data_qc_gDNA_long_qc$cdr3_seq_aa_q,sep="")

###To obtain the clonal persistant 
clone_df<-data.frame(table(data_qc_gDNA_long_qc$igh_clone_id,data_qc_gDNA_long_qc$specimen_label))
colnames(clone_df)<-c("clone","specimen","count")
clone_df$clin = reads_clones_annot_Long_qc[clone_df$specimen,1]
clone_df$time = reads_clones_annot_Long_qc[clone_df$specimen,4]
clone_df$sample = reads_clones_annot_Long_qc[clone_df$specimen,2]

clone_df_noceros = clone_df[which(clone_df$count!=0),]


g1<-ggplot(clone_df[which(clone_df$clin=="NP"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clonal persistant")

g2<-ggplot(clone_df_noceros[which(clone_df_noceros$clin=="PNR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clonal persistant")

g3<-ggplot(clone_df_noceros[which(clone_df_noceros$clin=="PR"),], aes(x=time,y=count,group=clone,fill=clone)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=clone)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "Clonal persistant")

tiff("Clonal_persistant_gDNA_long.tiff",res=300,h=1700,w=2700)
multiplot(g1,g2,g3)
dev.off()


##### V gene usage plots ####

############
### gDNA ###
############
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
#data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

### longitudinal
data_qc_gDNA_long<-data_qc_gDNA[which(data_qc_gDNA$clin!="AR" & data_qc_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_gDNA>100),]

id<-match(data_qc_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_gDNA_long_qc<-data_qc_gDNA_long[which(is.na(id)==F),]
data_qc_gDNA_long_qc$specimen_label<-factor(data_qc_gDNA_long_qc$specimen_label)

v_gene_df<-data.frame(table(data_qc_gDNA_long_qc$v_gene,data_qc_gDNA_long_qc$specimen_label))
colnames(v_gene_df)<-c("v_gene","specimen","count")
v_gene_df$clin = reads_clones_annot_Long_qc[v_gene_df$specimen,1] 
v_gene_df$time = reads_clones_annot_Long_qc[v_gene_df$specimen,4] 
v_gene_df$sample = reads_clones_annot_Long_qc[v_gene_df$specimen,2]
v_gene_df$clones = reads_clones_annot_Long_qc[v_gene_df$specimen,43]
v_gene_df$freq<-v_gene_df$count/v_gene_df$clones

g1<-ggplot(v_gene_df[which(v_gene_df$clin=="NP"),], aes(x=time,y=count,group=v_gene,fill=v_gene)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  facet_grid(clin ~ sample) + labs(x = "time", y = "V gene usage")

g2<-ggplot(v_gene_df[which(v_gene_df$clin=="PNR"),], aes(x=time,y=freq,group=v_gene,fill=v_gene)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  scale_y_continuous(limit = c(0,1)) + facet_grid(clin ~ sample) + labs(x = "time", y = "V gene usage")

g3<-ggplot(v_gene_df[which(v_gene_df$clin=="PR"),], aes(x=time,y=freq,group=v_gene,fill=v_gene)) +
  scale_x_continuous(breaks = c(0,6,12,24,32)) + geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  scale_y_continuous(limit = c(0,1)) + facet_grid(clin ~ sample) + labs(x = "time", y = "V gene usage")

tiff("v_gene_usage_gDNA_long.tiff",res=300,h=1700,w=2700)
multiplot(g1,g2,g3)
dev.off()


### AR
data_qc_gDNA_AR<-data_qc_gDNA[which(data_qc_gDNA$clin=="AR" | data_qc_gDNA$clin=="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]
reads_clones_annot_AR_qc<-reads_clones_annot_AR[which(reads_clones_annot_AR$reads_gDNA>100),]

id<-match(data_qc_gDNA_AR$specimen_label,reads_clones_annot_AR_qc$specimen_id)
data_qc_gDNA_AR_qc<-data_qc_gDNA_AR[which(is.na(id)==F),]
data_qc_gDNA_AR_qc$specimen_label<-factor(data_qc_gDNA_AR_qc$specimen_label)

v_gene_df<-data.frame(table(data_qc_gDNA_AR_qc$v_gene,data_qc_gDNA_AR_qc$specimen_label))
colnames(v_gene_df)<-c("v_gene","specimen","count")
v_gene_df$clin = reads_clones_annot_AR_qc[v_gene_df$specimen,1] 
v_gene_df$time = reads_clones_annot_AR_qc[v_gene_df$specimen,4] 
v_gene_df$sample = reads_clones_annot_AR_qc[v_gene_df$specimen,2]
v_gene_df$clones = reads_clones_annot_AR_qc[v_gene_df$specimen,43]
v_gene_df$reads = reads_clones_annot_AR_qc[v_gene_df$specimen,35]
v_gene_df$clin<-relevel(v_gene_df$clin,ref="pre-AR")
xx<-aggregate(v_gene_df$count, by=list(v_gene=v_gene_df$v_gene), FUN=sum)
v_gene_df$total = xx[v_gene_df$v_gene,2]
v_gene_df$freq<-v_gene_df$count/v_gene_df$total

tiff("v_gene_usage_gDNA_AR.tiff",res=300,h=1700,w=2700)
ggplot(v_gene_df, aes(x=clin,y=freq,group=v_gene,fill=v_gene)) +
  geom_area(aes(fill=v_gene)) + theme(legend.position="none") + 
  facet_grid(~ sample) + labs( y = "V gene usage")
dev.off()


############
### cDNA ###
############
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

### longitudinal
data_qc_cDNA_long<-data_qc_cDNA[which(data_qc_cDNA$clin!="AR" & data_qc_cDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_cDNA>100),]

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
  write.table(data_qc_cDNA_long_qc_specimen,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/data_qc_cDNA_long_qc_",i,".txt",sep=""),sep="\t",row.names = F)
  
  #assign(paste("vertex",i,sep=""),data.frame(table(data_qc_cDNA_long_qc_specimen$clone_CDR3)))
  #edges<-data.frame(table(data_qc_cDNA_long_qc_specimen$clone_CDR3,data_qc_cDNA_long_qc_specimen$igh_clone_id))
  #edges<-edges[which(edges$Freq==1),]
  #xx<-which(table(edges$Var2)>1)
  #id<-match(edges$Var2,names(xx))
  #edges<-edges[which(is.na(id)==F),]
  #write.table(edges,paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/edges",i,".txt",sep=""),sep="\t",row.names = F)
}

###Once is in network format I plot the network

##NP
files <- list.files("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/NP/")
for(i in files) {
  cat(i, "\n")
  edges <- read.delim(paste("/Users/Pinedasans/VDJ/Data/reads_by_sample/network_long_cDNA/NP/",i, sep = ""))
  j<-substr(i,6,10)
  nodes<-get(paste("vertex",j,sep=""))
  net<-graph_from_data_frame(d=edges,vertices = nodes,directed=T)
  V(net)$size <- V(net)$Freq/4
  V(net)$color <- c("chartreuse4")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  tiff(paste("/Users/Pinedasans/VDJ/Results/Network/network_",j,".tiff",sep=""),res=300,h=1500,w=1500)
  plot(net,vertex.label=NA,edge.curved=0.1,layout=layout_with_graphopt(net))
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




