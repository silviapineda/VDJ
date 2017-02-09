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
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: November, 2016
############################################################################################
library(gplots)
require(graphics)
library(RColorBrewer)
library("bcRep")

setwd("/Users/Pinedasans/Data/VDJ/")

##Read all the files and save into and Rdata all together
# files <- list.files("/Users/Pinedasans/Data/VDJ/part_tables_together/")
# 
# data <- c()
# for(i in files) {
#   cat(i, "\n")
#   t <- read.delim(paste("part_tables_together/", i, sep = ""))
#   data <- rbind(data, t)
# }
#save(data, file="new_merged_data.RData")

load("/Users/Pinedasans/Data/VDJ/new_merged_data.RData")

##Read the clinical annotations
clin_annot <- read.csv("/Users/Pinedasans/Data/VDJ/clin_annot_allsamples.csv")
rownames(clin_annot) <- clin_annot[,5]

##Some sample statistics
table(clin_annot$clin[!duplicated(clin_annot$Sample_id)])
unique(data[,c("specimen_label","amplification_label", "sample_label","specimen_description","participant_label")])

##############################
###Some Quality Control
#############################

##Discard all the V_score < 140
data_qc<-data[which(data$v_score>=140),]
##Discard the non-functional sequences
data_qc<-data_qc[which(data_qc$productive=="t"),]

#############################################
####Prepare data without downsampling#######
############################################

data_qc$clin = clin_annot[data_qc$specimen_label,1] ###Add the type of clinical endpoint
data_qc$time = clin_annot[data_qc$specimen_label,4] ###Add the time it was taking
data_qc$sample_id = clin_annot[data_qc$specimen_label,2] ###Add the sample_id

##Extract the gene from the segment with the allele
data_qc$v_gene <- gsub("\\*", "", substr(data_qc$v_segment, 1, 8))
data_qc$j_gene <- gsub("\\*", "", substr(data_qc$j_segment, 1, 5))
data_qc$d_gene <- gsub("\\*", "", substr(data_qc$d_segment, 1, 8))

###Extrac the CDR3 region
data_qc$cdr3_seq <- gsub(" ","", data_qc$cdr3_seq_nt_q)
###Extract the isotype
data_qc$isotype <- substr(data_qc$isosubtype, 1, 4)

##Read counts and clones per sample and data point
read_count <- table(data_qc$specimen_label)
read_count_amplification <- table(data_qc$specimen_label,data_qc$amplification_template)
read_count_isotype <- table(data_qc$specimen_label, data_qc$isotype)
colnames(read_count_isotype)[1] = "UNMAPPED"
colnames(read_count_isotype) <- paste(colnames(read_count_isotype), "isotypes", sep = "_")

reads <- cbind(read_count,read_count_amplification,read_count_isotype)
colnames(reads)[1:6] <- c("total_reads","cDNA_reads","gDNA_reads")
 

##Count number of clones per sample using Krishna estimates
data_qc$V_J_lenghCDR3_Clone_igh = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),data_qc$igh_clone_id,sep="_")
read_count_ighClones<- unique(data_qc[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template")])
clones_igh<-data.matrix(table(read_count_ighClones$specimen_label,read_count_ighClones$amplification_template))
colnames(clones_igh)<-c("clones_cDNA","clones_gDNA")

reads_clones_igh<-cbind(reads,clones_igh)

###To obtaion the overlapping samples between clinical annotation and data
id.sample <- match(rownames(reads_clones_igh),clin_annot$specimen_id)
reads_clones_annot <- cbind(clin_annot[id.sample,], reads_clones_igh)
write.csv(reads_clones_annot, "total_reads_clones.csv", row.names = F)

save(data_qc,reads_clones_annot,file="/Users/Pinedasans/Data/VDJ/VDJ.Rdata")

#############################################
####Prepare data with downsampling#######
############################################

#####################################################
###  Identification of clones with python script ####
#####################################################
###Count the number of clones
#Clonal inference works by requiring the same V and J segments (not including the allele), 
#same CDR3 length, 
#and 90% nucleotide identity between CDR3s.

##To obtain the number of clones using the program I made in python to obtain the clonal Inference using downsampling by individuals
data_qc$V_J_lenghCDR3 = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),sep="_")
data_clonesInference<-data_qc[,c("specimen_label","sample_id","amplification_template","isosubtype","v_segment","d_segment",
                                 "j_segment","igh_clone_id","cdr3_seq_aa_q","clin","time","v_gene","j_gene","d_gene","V_J_lenghCDR3","cdr3_seq")]
data_clonesInference<-data_clonesInference[which(nchar(data_qc$cdr3_seq)!=0),]
write.table(data_clonesInference,"/Users/Pinedasans/Data/VDJ/data_clonesInference.txt",row.names = F,sep="\t") 

#Read data from the output file from CloneInference.py
##Read all the files and save into and Rdata all together
files <- list.files("/Users/Pinedasans/Data/VDJ/ClonesInfered/")
data_clonesInference <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.csv(paste("ClonesInfered/", i, sep = ""))
  t$dataset<-rep(i,nrow(t))
  data_clonesInference <- rbind(data_clonesInference, t)
}
save(data_clonesInference, file="ClonesInfered_downsampled.RData")

##Read counts and clones per sample and data point
read_count <- table(data_clonesInference$specimen_label,data_clonesInference$dataset)
read_count<-apply(read_count,1,mean)

read_count_gDNA<-NULL
read_count_cDNA<-NULL
for (i in 1:5){
  read_count_gDNA<- rbind(read_count_gDNA,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                    data_clonesInference$amplification_template[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  read_count_cDNA<- rbind(read_count_cDNA,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference$amplification_template[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
}
read_count_gDNA<-apply(read_count_gDNA,2,mean)
read_count_cDNA<-apply(read_count_cDNA,2,mean)
read_count_amplification<-cbind(read_count_gDNA,read_count_cDNA)

data_clonesInference$isotype <- substr(data_clonesInference$isosubtype, 1, 4)
read_count_UNMAPPED<-NULL
read_count_IGHA<-NULL
read_count_IGHD<-NULL
read_count_IGHE<-NULL
read_count_IGHG<-NULL
read_count_IGHM<-NULL
for (i in 1:5){
  read_count_UNMAPPED<- rbind(read_count_UNMAPPED,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference$isotype[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
  read_count_IGHA<- rbind(read_count_IGHA,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                    data_clonesInference$isotype[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  read_count_IGHD<- rbind(read_count_IGHD,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference$isotype[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,3])
  read_count_IGHE<- rbind(read_count_IGHE,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference$isotype[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,4])
  read_count_IGHG<- rbind(read_count_IGHG,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference$isotype[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,5])
  read_count_IGHM<- rbind(read_count_IGHM,table(data_clonesInference$specimen_label[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference$isotype[which(data_clonesInference$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,6])
  
}
read_count_UNMAPPED<-apply(read_count_UNMAPPED,2,mean)
read_count_IGHA<-apply(read_count_IGHA,2,mean)
read_count_IGHD<-apply(read_count_IGHD,2,mean)
read_count_IGHE<-apply(read_count_IGHE,2,mean)
read_count_IGHG<-apply(read_count_IGHG,2,mean)
read_count_IGHM<-apply(read_count_IGHM,2,mean)

read_count_isotype<-cbind(read_count_UNMAPPED,read_count_IGHA,read_count_IGHD,read_count_IGHE,read_count_IGHG,read_count_IGHM)
colnames(read_count_isotype)[1] = "UNMAPPED"
colnames(read_count_isotype) <- paste(colnames(read_count_isotype), "isotypes", sep = "_")

reads <- cbind(read_count,read_count_amplification,read_count_isotype)

##Count number of clones per sample downsampled 
data_clonesInference$V_J_lenghCDR3_Clone = paste(data_clonesInference$V_J_lenghCDR3,data_clonesInference$numberClone,sep="_")
read_count_Clones<- unique(data_clonesInference[,c("specimen_label","V_J_lenghCDR3_Clone","amplification_template","dataset")])
clones_down<-data.matrix(table(read_count_Clones$specimen_label,read_count_Clones$dataset))
clones_down<-apply(clones_down,1,mean)

clones_down_gDNA<-NULL
clones_down_cDNA<-NULL
for (i in 1:5){
  clones_down_gDNA<- rbind(clones_down_gDNA,table(read_count_Clones$specimen_label[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                  read_count_Clones$amplification_template[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  clones_down_cDNA<- rbind(clones_down_cDNA,table(read_count_Clones$specimen_label[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                  read_count_Clones$amplification_template[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
}
clones_down_gDNA<-apply(clones_down_gDNA,2,mean)
clones_down_cDNA<-apply(clones_down_cDNA,2,mean)
clones_down_amplification<-cbind(clones_down_gDNA,clones_down_cDNA)
clones_down<-cbind(clones_down,clones_down_amplification)

####Cpunt number of cloned per sample igh
data_clonesInference$V_J_lenghCDR3_Clone_igh = paste(data_clonesInference$V_J_lenghCDR3,data_clonesInference$igh_clone_id,sep="_")
read_count_Clones_igh<- unique(data_clonesInference[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template","dataset")])
clones_down_igh<-data.matrix(table(read_count_Clones_igh$specimen_label,read_count_Clones_igh$dataset))
clones_down_igh<-apply(clones_down_igh,1,mean)

clones_down_igh_gDNA<-NULL
clones_down_igh_cDNA<-NULL
for (i in 1:5){
  clones_down_igh_gDNA<- rbind(clones_down_igh_gDNA,table(read_count_Clones_igh$specimen_label[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                          read_count_Clones_igh$amplification_template[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  clones_down_igh_cDNA<- rbind(clones_down_igh_cDNA,table(read_count_Clones_igh$specimen_label[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                          read_count_Clones_igh$amplification_template[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
}
clones_down_igh_gDNA<-apply(clones_down_igh_gDNA,2,mean)
clones_down_igh_cDNA<-apply(clones_down_igh_cDNA,2,mean)
clones_down_igh_amplification<-cbind(clones_down_igh_gDNA,clones_down_igh_cDNA)

clones_down_igh<-cbind(clones_down_igh,clones_down_igh_amplification)

reads_clones_down <- cbind(reads,clones_down,clones_down_igh)

###To obtaion the overlapping samples between clinical annotation and data
id.sample <- match(rownames(reads_clones_down),clin_annot$specimen_id)
reads_clones_down_annot <- cbind(clin_annot[id.sample,], reads_clones_down)
write.csv(reads_clones_down_annot, "total_reads_clones_down.csv", row.names = F)

save(data_clonesInference,reads_clones_down_annot,file="/Users/Pinedasans/Data/VDJ/VDJ_downSampled.Rdata")



#######################################################################
######## Downsampling the data directly without re-doing clones  ######
######################################################################
#Downsampling  

##gDNA 1062 reads is the minimum (after QC) for LongData
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
specimen_unique<-unique(data_qc_gDNA$specimen_label)

clones_igh_gDNA<-matrix(NA,length(specimen_unique),100)
for (j in 1:100){
  for (i in 1:length(specimen_unique)){
    data_specimen_unique<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(nrow(data_specimen_unique)>=1062){
      data_specimen_unique.down<-data_specimen_unique[sample(nrow(data_specimen_unique),1062),]
      clones_igh_gDNA[i,j]<- length(unique(data_specimen_unique.down[,"V_J_lenghCDR3_Clone"]))
    }
  }
}

clones_igh_gDNA_mean<-as.integer(apply(clones_igh_gDNA,1,mean))
names(clones_igh_gDNA_mean)<-specimen_unique
id.specimen<-match(reads_clones_annot$specimen_id,names(clones_igh_gDNA_mean))
reads_clones_annot$clones_igh_gDNA_downsample<-clones_igh_gDNA_mean[id.specimen]

##gDNA 5124 reads is the minimum (after QC) for ARdata
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
specimen_unique<-unique(data_qc_gDNA$specimen_label)

clones_igh_gDNA<-matrix(NA,length(specimen_unique),100)
for (j in 1:100){
  for (i in 1:length(specimen_unique)){
    data_specimen_unique<-data_qc_gDNA[which(data_qc_gDNA$specimen_label==specimen_unique[i]),]
    if(nrow(data_specimen_unique)>=5124){
      data_specimen_unique.down<-data_specimen_unique[sample(nrow(data_specimen_unique),5124),]
      clones_igh_gDNA[i,j]<- length(unique(data_specimen_unique.down[,"V_J_lenghCDR3_Clone"]))
    }
  }
}

clones_igh_gDNA_mean<-as.integer(apply(clones_igh_gDNA,1,mean))
names(clones_igh_gDNA_mean)<-specimen_unique
id.specimen<-match(reads_clones_annot$specimen_id,names(clones_igh_gDNA_mean))
reads_clones_annot$clones_igh_gDNA_downsample_AR<-clones_igh_gDNA_mean[id.specimen]

##cDNA 62173 read is the minimum (after QC)
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
specimen_unique<-unique(data_qc_cDNA$specimen_label)

clones_igh_cDNA<-matrix(NA,length(specimen_unique),100)
for (j in 1:100){
  for (i in 1:length(specimen_unique)){
    data_specimen_unique<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen_unique[i]),]
    if(nrow(data_specimen_unique)>=62173){
      data_specimen_unique.down<-data_specimen_unique[sample(nrow(data_specimen_unique),62173),]
      clones_igh_cDNA[i,j]<- length(unique(data_specimen_unique.down[,"V_J_lenghCDR3_Clone"]))
    }
  }
}

clones_igh_cDNA_mean<-as.integer(apply(clones_igh_cDNA,1,mean))
names(clones_igh_cDNA_mean)<-specimen_unique
id.specimen<-match(reads_clones_annot$specimen_id,names(clones_igh_cDNA_mean))
reads_clones_annot$clones_igh_cDNA_downsample<-clones_igh_cDNA_mean[id.specimen]

