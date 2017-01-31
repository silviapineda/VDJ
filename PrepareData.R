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
clin_annot <- read.csv("/Users/Pinedasans/Data/VDJ/clin_annot.csv")
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

####Prepare data
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

##Read counts per sample and data point
read_count <- table(data_qc$specimen_label)
read_count_amplification <- table(data_qc$specimen_label,data_qc$amplification_template)
read_count_isotype <- table(data_qc$specimen_label, data_qc$isotype)
read_count_run <- table(data_qc$specimen_label,data_qc$run_label)
colnames(read_count_isotype)[1] = "UNMAPPED"
colnames(read_count_isotype) <- paste(colnames(read_count_isotype), "isotypes", sep = "_")

reads <- cbind(read_count,read_count_amplification,read_count_run,read_count_isotype)
colnames(reads)[1:6] <- c("total_reads","cDNA_reads","gDNA_reads","M154_reads","M155_reads","T7_reads")
  

###Identification of clones
###Count the number of clones
#Clonal inference works by requiring the same V and J segments (not including the allele), 
#same CDR3 length, 
#and 90% nucleotide identity between CDR3s. 
data_qc$V_J_lenghCDR3 = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),sep="_")
data_qc$V_J_lenghCDR3_Clone = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),data_qc$igh_clone_id,sep="_")

##Count number of clones per sample using Krishna estimates
read_count_ighClones<- unique(data_qc[,c("specimen_label","V_J_lenghCDR3_Clone")])
clones_igh<-data.matrix(table(read_count_ighClones$specimen_label))

##To obtain the number of clones using the program I made in python to obtain the clonal Inference
data_clonesInference<-data_qc[,c("specimen_label","V_J_lenghCDR3","cdr3_seq")]
data_clonesInference<-data_clonesInference[which(nchar(data_qc$cdr3_seq)!=0),]
write.table(data_clonesInference,"/Users/Pinedasans/Data/VDJ/data_clonesInference.txt",row.names = F,sep="\t") 

###Considering only cDNA
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
data_clonesInference_cDNA<-data_qc_cDNA[,c("specimen_label","V_J_lenghCDR3","cdr3_seq")]
data_clonesInference_cDNA<-data_clonesInference_cDNA[which(nchar(data_qc_cDNA$cdr3_seq)!=0),]
write.table(data_clonesInference_cDNA,"/Users/Pinedasans/Data/VDJ/data_clonesInference_cDNA.txt",row.names = F,sep="\t") 

###Considering only gDNA
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
data_clonesInference_gDNA<-data_qc_gDNA[,c("specimen_label","V_J_lenghCDR3","cdr3_seq")]
data_clonesInference_gDNA<-data_clonesInference_gDNA[which(nchar(data_qc_gDNA$cdr3_seq)!=0),]
write.table(data_clonesInference_gDNA,"/Users/Pinedasans/Data/VDJ/data_clonesInference_gDNA.txt",row.names = F,sep="\t") 

#Clone count per sample
#Read data from the output file from CloneInference.py
clone_count <- read.csv("/Users/Pinedasans/Data/VDJ/clones_count_sample.csv")
clone_count <- clone_count[,2:4]
id.sample <- match(rownames(reads),clone_count$specimen_label)
reads_clones <- cbind(reads,clones_igh[,1],clone_count[id.sample,1:2])
colnames(reads_clones)[13:14]<-c("clones_igh","total_clones")

###To obtaion the overlapping samples between clinical annotation and data
id.sample <- match(rownames(reads_clones),clin_annot$specimen_id)

reads_clones_annot <- cbind(clin_annot[id.sample,], reads_clones)
write.csv(reads_clones_annot, "total_reads_clones_new.csv", row.names = F)

#####Only using cDNA reads
clone_count_cDNA <- read.csv("/Users/Pinedasans/Data/VDJ/clones_count_sample_cDNA.csv")
clone_count_cDNA <- clone_count_cDNA[,2:4]
id.sample <- match(rownames(reads),clone_count_cDNA$specimen_label)
reads_clones_cDNA <- cbind(reads,clones_igh[,1],clone_count_cDNA[id.sample,1:2])
colnames(reads_clones_cDNA)[13:14]<-c("clones_igh","total_clones")
id.sample <- match(rownames(reads_clones_cDNA),clin_annot$specimen_id)
reads_clones_annot_cDNA <- cbind(clin_annot[id.sample,], reads_clones_cDNA)
write.csv(reads_clones_annot_cDNA, "total_reads_clones_new_cDNA.csv", row.names = F)

#####Only using gDNA reads
clone_count_gDNA <- read.csv("/Users/Pinedasans/Data/VDJ/clones_count_sample_gDNA.csv")
clone_count_gDNA <- clone_count_gDNA[,2:4]
id.sample <- match(rownames(reads),clone_count_gDNA$specimen_label)
reads_clones_gDNA <- cbind(reads,clones_igh[,1],clone_count_gDNA[id.sample,1:2])
colnames(reads_clones_gDNA)[13:14]<-c("clones_igh","total_clones")
id.sample <- match(rownames(reads_clones_gDNA),clin_annot$specimen_id)
reads_clones_annot_gDNA <- cbind(clin_annot[id.sample,], reads_clones_gDNA)
write.csv(reads_clones_annot_gDNA, "total_reads_clones_new_gDNA.csv", row.names = F)


save(data_qc,reads_clones_annot,reads_clones_annot_cDNA,reads_clones_annot_gDNA,file="/Users/Pinedasans/Data/VDJ/VDJ.Rdata")
