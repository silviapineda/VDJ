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

setwd("/Users/Pinedasans/VDJ/Data/")

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

###This is the data received from Krishna
load("/Users/Pinedasans/VDJ/Data/new_merged_data.RData")

##Read the clinical annotations
clin_annot <- read.csv("/Users/Pinedasans/VDJ/Data/clin_annot_allsamples.csv")
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
colnames(reads)[1:3] <- c("total_reads","reads_cDNA","reads_gDNA")
 

##Count number of clones per sample using Krishna estimates
data_qc$V_J_lenghCDR3_Clone_igh = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),data_qc$igh_clone_id,sep="_")
read_count_ighClones<- unique(data_qc[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template")])
clones_igh<-data.matrix(table(read_count_ighClones$specimen_label,read_count_ighClones$amplification_template))
colnames(clones_igh)<-c("clones_cDNA","clones_gDNA")

reads_clones_igh<-cbind(reads,clones_igh)

###To obtain the overlapping samples between clinical annotation and data
id.sample <- match(rownames(reads_clones_igh),clin_annot$specimen_id)
reads_clones_annot <- cbind(clin_annot[id.sample,], reads_clones_igh)
write.csv(reads_clones_annot, "/Users/Pinedasans/VDJ/Data/total_reads_clones.csv", row.names = F)

save(data_qc,reads_clones_annot,file="/Users/Pinedasans/VDJ/Data//VDJ.Rdata")




