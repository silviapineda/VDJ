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
files <- list.files("/Users/Pinedasans/VDJ/Data/part_tables_together")

data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("part_tables_together/", i, sep = ""))
  data <- rbind(data, t)
}
save(data, file="new_merged_data.RData")

###This is the data received from Krishna
load("/Users/Pinedasans/VDJ/Data/new_merged_data.RData")

##plot histogram for the V score
hist(data$v_score,xlab = "V gene score", main="Histogram: Scott Data")
abline(v=140,col="red")

##Read the clinical annotations
clin_annot <- read.csv("/Users/Pinedasans/VDJ/Data/clin_annot_allsamples.csv")
rownames(clin_annot) <- clin_annot[,5]

##Some sample statistics
table(clin_annot$clin[!duplicated(clin_annot$Sample_id)])
unique(data[,c("specimen_label","amplification_label", "sample_label","specimen_description","participant_label")])

##############################
###Some Quality Control
#############################
##Discard the non-functional sequences
data_qc<-data[which(data$productive=="t"),]

##Discard all the V_score < 140
data_qc<-data[which(data$v_score>=140),]
##Discard the non-functional sequences
data_qc<-data_qc[which(data_qc$productive=="t"),]

#########################
####Prepare data #######
#######################
data_qc_order = data_qc[order(data_qc$specimen_label),]
data_qc_order$clin = clin_annot[data_qc_order$specimen_label,1] ###Add the type of clinical endpoint
data_qc_order$time = clin_annot[data_qc_order$specimen_label,4] ###Add the time it was taking
data_qc_order$sample_id = clin_annot[data_qc_order$specimen_label,2] ###Add the sample_id

##Extract the gene from the segment with the allele
data_qc_order$v_gene <- gsub("\\*", "", substr(data_qc_order$v_segment, 1, 8))
data_qc_order$j_gene <- gsub("\\*", "", substr(data_qc_order$j_segment, 1, 5))
data_qc_order$d_gene <- gsub("\\*", "", substr(data_qc_order$d_segment, 1, 8))

###Extract the CDR3 region
data_qc_order$cdr3_seq <- gsub(" ","", data_qc_order$cdr3_seq_nt_q)
###Extract the isotype
data_qc_order$isotype <- substr(data_qc_order$isosubtype, 1, 4)


##Count the somatic hypermutations
data_qc_order$Vlength<-nchar(as.character(data_qc_order$v_sequence))
v_sequence = as.character(data_qc_order$v_sequence)
data_qc_order$SHM<-sapply(regmatches(v_sequence, gregexpr("[A-Z]", v_sequence, perl=TRUE)), length)
data_qc_order$SHM_freq<-data_qc_order$SHM/data_qc_order$Vlength

##count the CDR3 length
data_qc_order$CDR3_length<-nchar(as.character(data_qc_order$cdr3_seq))

##count the read length
data_qc_order$read_length<-nchar(as.character(data_qc_order$trimmed_sequence))

##Read counts and clones per sample and data point
read_count <- table(data_qc_order$specimen_label)
read_count_amplification <- table(data_qc_order$specimen_label,data_qc_order$amplification_template)
read_count_isotype <- table(data_qc_order$specimen_label, data_qc_order$isotype)
colnames(read_count_isotype)[1] = "UNMAPPED"
colnames(read_count_isotype) <- paste(colnames(read_count_isotype), "isotypes", sep = "_")

reads <- cbind(read_count,read_count_amplification,read_count_isotype)
colnames(reads)[1:3] <- c("total_reads","reads_cDNA","reads_gDNA")
 

##Count number of clones per sample using Krishna estimates
data_qc_order$V_J_lenghCDR3_Clone_igh = paste(data_qc_order$v_gene, data_qc_order$j_gene, nchar(data_qc_order$cdr3_seq),data_qc_order$igh_clone_id,sep="_")
read_count_ighClones<- unique(data_qc_order[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template")])
clones_igh<-data.matrix(table(read_count_ighClones$specimen_label,read_count_ighClones$amplification_template))
colnames(clones_igh)<-c("clones_cDNA","clones_gDNA")

reads_clones_igh<-cbind(reads,clones_igh)

#####SHM
sum_SHM_freq<-aggregate(data_qc_order$SHM_freq, by=list(data_qc_order$specimen_label,data_qc_order$amplification_template), FUN=sum)
sum_SHM_gDNA<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="gDNA"),]
sum_SHM_cDNA<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="cDNA"),]
id_gDNA<-match(rownames(reads_clones_igh),sum_SHM_gDNA$Group.1)
id_cDNA<-match(rownames(reads_clones_igh),sum_SHM_cDNA$Group.1)
reads_clones_igh_SHM<-cbind(reads_clones_igh,sum_SHM_gDNA$x[id_gDNA],sum_SHM_cDNA$x[id_cDNA])
colnames(reads_clones_igh_SHM)[12:13]<-c("SHM_gDNA","SHM_cDNA")

###length of CDR3
cdr3_length<-aggregate(data_qc_order$CDR3_length,by=list(data_qc_order$specimen_label,data_qc_order$amplification_template), FUN=mean)
cdr3_length_gDNA<-cdr3_length[which(cdr3_length$Group.2=="gDNA"),]
cdr3_length_cDNA<-cdr3_length[which(cdr3_length$Group.2=="cDNA"),]
id_gDNA<-match(rownames(reads_clones_igh),cdr3_length_gDNA$Group.1)
id_cDNA<-match(rownames(reads_clones_igh),cdr3_length_cDNA$Group.1)
reads_clones_igh_cdr3length<-cbind(reads_clones_igh_SHM,cdr3_length_gDNA$x[id_gDNA],cdr3_length_cDNA$x[id_cDNA])
colnames(reads_clones_igh_cdr3length)[14:15]<-c("mean_CDR3_length_gDNA","mean_CDR3_length_cDNA")


###To obtain the overlapping samples between clinical annotation and data
id.sample <- match(rownames(reads_clones_igh_cdr3length),clin_annot$specimen_id)
reads_clones_annot <- cbind(clin_annot[id.sample,], reads_clones_igh_cdr3length)
#write.csv(reads_clones_annot, "/Users/Pinedasans/VDJ/Data/total_reads_clones.csv", row.names = F)
data_qc<-data_qc_order
write.csv(reads_clones_annot,"reads_clones_annotation.csv")

##Read the annotation file with the new subject_id
reads_clones_annot<-read.csv("/Users/Pinedasans/VDJ/Data/reads_clones_annotation.csv")

save(data_qc,reads_clones_annot,file="/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")

###save the data to call the clones by all samples
data_qc$unique_id<-seq(1,nrow(data_qc))
data_clonesInference<-data_qc[,c("unique_id","specimen_label","sample_id","amplification_template","isosubtype","v_segment","d_segment",
                                 "j_segment","trimmed_sequence","replicate_label","v_sequence","d_sequence","j_sequence","igh_clone_id","cdr3_seq_aa_q","clin","time","v_gene","j_gene","d_gene","isotype",
                                 "Vlength", "SHM",  "SHM_freq", "CDR3_length","V_J_lenghCDR3","cdr3_seq")]
data_clonesInference<-data_clonesInference[which(nchar(data_qc$cdr3_seq)!=0),]
write.table(data_clonesInference,file="/Users/Pinedasans/VDJ/Data/data_clonesInference.txt",row.names = F,sep="\t")


###Merge nucleotides output with original data
load(file="/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")
data_qc<-data_qc[which(nchar(data_qc$cdr3_seq)!=0),]
nucleotides_output<-read.csv("/Users/Pinedasans/VDJ/Data/ClonesInferedAll_90.csv")
nucleotides_clones_unique<-unique(nucleotides_output[,c("V_J_lenghCDR3","cdr3_seq","CloneId")])
data_merge<-merge(data_qc,nucleotides_clones_unique,by=c("V_J_lenghCDR3","cdr3_seq"))
save(data_merge,reads_clones_annot,file="/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")



#### prepare downsampling data in cDNA ####
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")
reads_clones_annot_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]
reads_clones_annot_cDNA<-reads_clones_annot[which(reads_clones_annot$reads_cDNA>100),]
minread<-min(reads_clones_annot_cDNA$reads_cDNA) #62173
data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]


specimen<-reads_clones_annot_cDNA$specimen_id
clones_igh<-matrix(NA,103,10)
for ( j in 1:10){
  data_qc_cDNA_downsampling<-NULL
  for(i in 1:length(specimen)){
    print(i)
    data_qc_cDNA_specimen<-data_qc_cDNA[which(data_qc_cDNA$specimen_label==specimen[i]),]
    data_qc_cDNA_downsampling<-rbind(data_qc_cDNA_downsampling,data_qc_cDNA_specimen[sample(minread),])
  }

  read_count_ighClones<- unique(data_qc_cDNA_downsampling[,c("specimen_label","V_J_lenghCDR3_Clone_igh")])
  clones_igh[,j]<-data.matrix(table(read_count_ighClones$specimen_label))
}
