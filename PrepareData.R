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
colnames(reads)[1:3] <- c("total_reads","cDNA_reads","gDNA_reads")
 

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

##The downsampling is done in the python script
##To obtain the number of clones using the program I made in python to obtain the clonal Inference using downsampling by individuals
data_qc$V_J_lenghCDR3 = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),sep="_")
data_clonesInference<-data_qc[,c("specimen_label","sample_id","amplification_template","isosubtype","v_segment","d_segment",
                                 "j_segment","igh_clone_id","cdr3_seq_aa_q","clin","time","v_gene","j_gene","d_gene","V_J_lenghCDR3","cdr3_seq")]
data_clonesInference<-data_clonesInference[which(nchar(data_qc$cdr3_seq)!=0),]
write.table(data_clonesInference,"/Users/Pinedasans/Data/VDJ/data_clonesInference.txt",row.names = F,sep="\t") 

####Downsampling by cDNA and gDNA and put together to pass python program
##cDNA 62173 read is the minimum (after QC)
data_clonesInference_cDNA<-data_clonesInference[which(data_clonesInference$amplification_template=="cDNA"),]
specimen_unique<-unique(data_clonesInference_cDNA$specimen_label)
data_specimen_unique_cDNA.down<-NULL
for (i in 1:length(specimen_unique)){
  data_specimen_unique<-data_clonesInference_cDNA[which(data_clonesInference_cDNA$specimen_label==specimen_unique[i]),]
    if(nrow(data_specimen_unique)>=62173){
      data_specimen_unique_cDNA.down<-rbind(data_specimen_unique_cDNA.down,data_specimen_unique[sample(nrow(data_specimen_unique),62173),])
    }
}

##gDNA 1062 read is the minimum (after QC)
data_clonesInference_gDNA<-data_clonesInference[which(data_clonesInference$amplification_template=="gDNA"),]
specimen_unique<-unique(data_clonesInference_gDNA$specimen_label)
data_specimen_unique_gDNA.down<-NULL
for (i in 1:length(specimen_unique)){
  data_specimen_unique<-data_clonesInference_gDNA[which(data_clonesInference_gDNA$specimen_label==specimen_unique[i]),]
  if(nrow(data_specimen_unique)>=1062){
    data_specimen_unique_gDNA.down<-rbind(data_specimen_unique_gDNA.down,data_specimen_unique[sample(nrow(data_specimen_unique),1062),])
  }
}

data_clonesInference_amplification_down <- rbind(data_specimen_unique_cDNA.down,data_specimen_unique_gDNA.down)
write.table(data_clonesInference_amplification_down,"/Users/Pinedasans/Data/VDJ/data_clonesInference_amplification_down.txt",row.names = F,sep="\t") 


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

#############################################################
#### ReadData DownSampling by Individual SecondMinimum ######
#############################################################

#Read data from the output file from CloneInference.py
##Read all the files and save into and Rdata all together
files <- list.files("/Users/Pinedasans/Data/VDJ/ClonesInferedIndividual_2minimum/")
data_clonesInference_2 <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.csv(paste("ClonesInferedIndividual_2minimum/", i, sep = ""))
  t$dataset<-rep(i,nrow(t))
  data_clonesInference_2 <- rbind(data_clonesInference_2, t)
}
save(data_clonesInference_2, file="ClonesInfered_2_downsampled.RData")

##Read counts and clones per sample and data point
read_count <- table(data_clonesInference_2$specimen_label,data_clonesInference_2$dataset)
read_count<-apply(read_count,1,mean)

read_count_gDNA<-NULL
read_count_cDNA<-NULL
for (i in 1:5){
  read_count_gDNA<- rbind(read_count_gDNA,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_2$amplification_template[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  read_count_cDNA<- rbind(read_count_cDNA,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_2$amplification_template[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
}
read_count_gDNA<-apply(read_count_gDNA,2,mean)
read_count_cDNA<-apply(read_count_cDNA,2,mean)
read_count_amplification<-cbind(read_count_gDNA,read_count_cDNA)

data_clonesInference_2$isotype <- substr(data_clonesInference_2$isosubtype, 1, 4)
read_count_UNMAPPED<-NULL
read_count_IGHA<-NULL
read_count_IGHD<-NULL
read_count_IGHE<-NULL
read_count_IGHG<-NULL
read_count_IGHM<-NULL
for (i in 1:5){
  read_count_UNMAPPED<- rbind(read_count_UNMAPPED,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                        data_clonesInference_2$isotype[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
  read_count_IGHA<- rbind(read_count_IGHA,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_2$isotype[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  read_count_IGHD<- rbind(read_count_IGHD,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_2$isotype[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,3])
  read_count_IGHE<- rbind(read_count_IGHE,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_2$isotype[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,4])
  read_count_IGHG<- rbind(read_count_IGHG,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_2$isotype[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,5])
  read_count_IGHM<- rbind(read_count_IGHM,table(data_clonesInference_2$specimen_label[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_2$isotype[which(data_clonesInference_2$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,6])
  
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
data_clonesInference_2$V_J_lenghCDR3_Clone = paste(data_clonesInference_2$V_J_lenghCDR3,data_clonesInference_2$numberClone,sep="_")
read_count_Clones<- unique(data_clonesInference_2[,c("specimen_label","V_J_lenghCDR3_Clone","amplification_template","dataset")])
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
data_clonesInference_2$V_J_lenghCDR3_Clone_igh = paste(data_clonesInference_2$V_J_lenghCDR3,data_clonesInference_2$igh_clone_id,sep="_")
read_count_Clones_igh<- unique(data_clonesInference_2[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template","dataset")])
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
reads_clones_down_annot_2 <- cbind(clin_annot[id.sample,], reads_clones_down)
write.csv(reads_clones_down_annot_2, "total_reads_clones_down_2.csv", row.names = F)

save(data_clonesInference_2,reads_clones_down_annot_2,file="/Users/Pinedasans/Data/VDJ/VDJ_downSampled_2.Rdata")

##################################################
#### ReadData DownSampling by Amplification ######
#################################################

#Read data from the output file from CloneInference.py
##Read all the files and save into and Rdata all together
files <- list.files("/Users/Pinedasans/Data/VDJ/ClonesInferedAmplification/")
data_clonesInference_amp <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.csv(paste("ClonesInferedAmplification/", i, sep = ""))
  t$dataset<-rep(i,nrow(t))
  data_clonesInference_amp <- rbind(data_clonesInference_amp, t)
}
save(data_clonesInference_amp, file="ClonesInferedAmplification_downsampled.RData")

##Read counts and clones per sample and data point
read_count <- table(data_clonesInference_amp$specimen_label,data_clonesInference_amp$dataset)
read_count<-apply(read_count,1,mean)

read_count_gDNA<-NULL
read_count_cDNA<-NULL
for (i in 1:5){
  read_count_gDNA<- rbind(read_count_gDNA,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                data_clonesInference_amp$amplification_template[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,2])
  read_count_cDNA<- rbind(read_count_cDNA,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                data_clonesInference_amp$amplification_template[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,1])
}
read_count_gDNA<-apply(read_count_gDNA,2,mean)
read_count_cDNA<-apply(read_count_cDNA,2,mean)
read_count_amplification<-cbind(read_count_gDNA,read_count_cDNA)

data_clonesInference_amp$isotype <- substr(data_clonesInference_amp$isosubtype, 1, 4)
read_count_UNMAPPED<-NULL
read_count_IGHA<-NULL
read_count_IGHD<-NULL
read_count_IGHE<-NULL
read_count_IGHG<-NULL
read_count_IGHM<-NULL
for (i in 1:5){
  read_count_UNMAPPED<- rbind(read_count_UNMAPPED,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                        data_clonesInference_amp$isotype[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,1])
  read_count_IGHA<- rbind(read_count_IGHA,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                data_clonesInference_amp$isotype[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,2])
  read_count_IGHD<- rbind(read_count_IGHD,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                data_clonesInference_amp$isotype[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,3])
  read_count_IGHE<- rbind(read_count_IGHE,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                data_clonesInference_amp$isotype[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,4])
  read_count_IGHG<- rbind(read_count_IGHG,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                data_clonesInference_amp$isotype[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,5])
  read_count_IGHM<- rbind(read_count_IGHM,table(data_clonesInference_amp$specimen_label[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                data_clonesInference_amp$isotype[which(data_clonesInference_amp$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,6])
  
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
data_clonesInference_amp$V_J_lenghCDR3_Clone = paste(data_clonesInference_amp$V_J_lenghCDR3,data_clonesInference_amp$numberClone,sep="_")
read_count_Clones<- unique(data_clonesInference_amp[,c("specimen_label","V_J_lenghCDR3_Clone","amplification_template","dataset")])
clones_down<-data.matrix(table(read_count_Clones$specimen_label,read_count_Clones$dataset))
clones_down<-apply(clones_down,1,mean)

clones_down_gDNA<-NULL
clones_down_cDNA<-NULL
for (i in 1:5){
  clones_down_gDNA<- rbind(clones_down_gDNA,table(read_count_Clones$specimen_label[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                  read_count_Clones$amplification_template[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,2])
  clones_down_cDNA<- rbind(clones_down_cDNA,table(read_count_Clones$specimen_label[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                  read_count_Clones$amplification_template[which(read_count_Clones$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,1])
}
clones_down_gDNA<-apply(clones_down_gDNA,2,mean)
clones_down_cDNA<-apply(clones_down_cDNA,2,mean)
clones_down_amplification<-cbind(clones_down_gDNA,clones_down_cDNA)
clones_down<-cbind(clones_down,clones_down_amplification)

####Cpunt number of cloned per sample igh
data_clonesInference_amp$V_J_lenghCDR3_Clone_igh = paste(data_clonesInference_amp$V_J_lenghCDR3,data_clonesInference_amp$igh_clone_id,sep="_")
read_count_Clones_igh<- unique(data_clonesInference_amp[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template","dataset")])
clones_down_igh<-data.matrix(table(read_count_Clones_igh$specimen_label,read_count_Clones_igh$dataset))
clones_down_igh<-apply(clones_down_igh,1,mean)

clones_down_igh_gDNA<-NULL
clones_down_igh_cDNA<-NULL
for (i in 1:5){
  clones_down_igh_gDNA<- rbind(clones_down_igh_gDNA,table(read_count_Clones_igh$specimen_label[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                          read_count_Clones_igh$amplification_template[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,2])
  clones_down_igh_cDNA<- rbind(clones_down_igh_cDNA,table(read_count_Clones_igh$specimen_label[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))],
                                                          read_count_Clones_igh$amplification_template[which(read_count_Clones_igh$dataset==paste("ClonesInfered_downsampled_byAmplification_",i,".csv",sep=""))])[,1])
}
clones_down_igh_gDNA<-apply(clones_down_igh_gDNA,2,mean)
clones_down_igh_cDNA<-apply(clones_down_igh_cDNA,2,mean)
clones_down_igh_amplification<-cbind(clones_down_igh_gDNA,clones_down_igh_cDNA)

clones_down_igh<-cbind(clones_down_igh,clones_down_igh_amplification)

reads_clones_down <- cbind(reads,clones_down,clones_down_igh)

###To obtaion the overlapping samples between clinical annotation and data
id.sample <- match(rownames(reads_clones_down),clin_annot$specimen_id)
reads_clones_down_annot_amp <- cbind(clin_annot[id.sample,], reads_clones_down)
write.csv(reads_clones_down_annot_amp, "total_reads_clones_down_amp.csv", row.names = F)

save(data_clonesInference_amp,reads_clones_down_annot_amp,file="/Users/Pinedasans/Data/VDJ/VDJ_amplification_downSampled.Rdata")

