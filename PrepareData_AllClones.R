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

###Load the merge data
load("/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")

##Read the clinical annotations
clin_annot <- read.csv("/Users/Pinedasans/VDJ/Data/clin_annot_allsamples.csv")
rownames(clin_annot) <- clin_annot[,5]


#########################
####Prepare data #######
#######################

##Read counts and clones per sample and data point
read_count <- table(data_merge$specimen_label)
read_count_amplification <- table(data_merge$specimen_label,data_merge$amplification_template)
read_count_isotype <- table(data_merge$specimen_label, data_merge$isotype)
colnames(read_count_isotype)[1] = "UNMAPPED"
colnames(read_count_isotype) <- paste(colnames(read_count_isotype), "isotypes", sep = "_")

data_merge$replicate_label<-as.character(data_merge$replicate_label)
data_merge$replicates_gDNA<-substr(data_merge$replicate_label, nchar(data_merge$replicate_label)-4, nchar(data_merge$replicate_label))
read_count_replicates <- table(data_merge$specimen_label, data_merge$replicates_gDNA)
read_count_replicates <- read_count_replicates[,1:6]

reads <- cbind(read_count,read_count_amplification,read_count_isotype,read_count_replicates)
colnames(reads)[1:3] <- c("total_reads","reads_cDNA","reads_gDNA")
 
##Count number of clones per sample using All samples
data_merge$V_J_lenghCDR3_CloneId = paste(data_merge$V_J_lenghCDR3,data_merge$CloneId,sep="_")
read_count_AllClones<- unique(data_merge[,c("specimen_label","V_J_lenghCDR3_CloneId","amplification_template")])
clones_All<-data.matrix(table(read_count_AllClones$specimen_label,read_count_AllClones$amplification_template))
colnames(clones_All)<-c("clones_cDNA","clones_gDNA")
##For isotypes
read_count_isotypes<- unique(data_merge[,c("specimen_label","V_J_lenghCDR3_CloneId","isotype")])
clones_isotypes<-data.matrix(table(read_count_isotypes$specimen_label,read_count_isotypes$isotype))
colnames(clones_isotypes)<-c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM")
##For replicates
read_count_replicates<- unique(data_merge[,c("specimen_label","V_J_lenghCDR3_CloneId","replicates_gDNA")])
clones_replicates<-data.matrix(table(read_count_replicates$specimen_label,read_count_replicates$replicates_gDNA))
clones_replicates<-clones_replicates[,1:6]
colnames(clones_replicates)<-c("clones_PCR_A","clones_PCR_B","clones_PCR_C","clones_PCR_D","clones_PCR_E","clones_PCR_F")

###Putting all together
reads_clones_all<-cbind(reads,clones_All,clones_isotypes,clones_replicates)

#####SHM
sum_SHM_freq<-aggregate(data_merge$SHM_freq, by=list(data_merge$specimen_label,data_merge$amplification_template), FUN=sum)
sum_SHM_gDNA<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="gDNA"),]
sum_SHM_cDNA<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="cDNA"),]
id_gDNA<-match(rownames(reads_clones_all),sum_SHM_gDNA$Group.1)
id_cDNA<-match(rownames(reads_clones_all),sum_SHM_cDNA$Group.1)
reads_clones_SHM<-cbind(sum_SHM_gDNA$x[id_gDNA],sum_SHM_cDNA$x[id_cDNA])
colnames(reads_clones_SHM)<-c("SHM_gDNA","SHM_cDNA")
##For isotypes
sum_SHM_freq<-aggregate(data_merge$SHM_freq, by=list(data_merge$specimen_label,data_merge$isotype), FUN=sum)
sum_SHM_unmapped<-sum_SHM_freq[which(sum_SHM_freq$Group.2==""),]
sum_SHM_IGHA<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHA"),]
sum_SHM_IGHD<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHD"),]
sum_SHM_IGHE<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHE"),]
sum_SHM_IGHG<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHG"),]
sum_SHM_IGHM<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHM"),]
sum_SHM_IGHM<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHM"),]

id_unmapped<-match(rownames(reads_clones_all),sum_SHM_unmapped$Group.1)
id_IGHA<-match(rownames(reads_clones_all),sum_SHM_IGHA$Group.1)
id_IGHD<-match(rownames(reads_clones_all),sum_SHM_IGHD$Group.1)
id_IGHE<-match(rownames(reads_clones_all),sum_SHM_IGHE$Group.1)
id_IGHG<-match(rownames(reads_clones_all),sum_SHM_IGHG$Group.1)
id_IGHM<-match(rownames(reads_clones_all),sum_SHM_IGHM$Group.1)

reads_clones_SHM_isotypes<-cbind(sum_SHM_unmapped$x[id_unmapped],sum_SHM_IGHA$x[id_IGHA],sum_SHM_IGHD$x[id_IGHD],
                        sum_SHM_IGHE$x[id_IGHE],sum_SHM_IGHG$x[id_IGHG],sum_SHM_IGHM$x[id_IGHM])
colnames(reads_clones_SHM_isotypes)<-c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHE","SHM_IGHG","SHM_IGHM")

##For repliactes
sum_SHM_freq<-aggregate(data_merge$SHM_freq, by=list(data_merge$specimen_label,data_merge$replicates_gDNA), FUN=sum)
sum_SHM_PCR_A<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="PCR_A"),]
sum_SHM_PCR_B<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="PCR_B"),]
sum_SHM_PCR_C<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="PCR_C"),]
sum_SHM_PCR_D<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="PCR_D"),]
sum_SHM_PCR_E<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="PCR_E"),]
sum_SHM_PCR_F<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="PCR_F"),]
id_PCR_A<-match(rownames(reads_clones_all),sum_SHM_PCR_A$Group.1)
id_PCR_B<-match(rownames(reads_clones_all),sum_SHM_PCR_B$Group.1)
id_PCR_C<-match(rownames(reads_clones_all),sum_SHM_PCR_C$Group.1)
id_PCR_D<-match(rownames(reads_clones_all),sum_SHM_PCR_D$Group.1)
id_PCR_E<-match(rownames(reads_clones_all),sum_SHM_PCR_E$Group.1)
id_PCR_F<-match(rownames(reads_clones_all),sum_SHM_PCR_F$Group.1)
reads_clones_SHM_replicates<-cbind(sum_SHM_PCR_A$x[id_PCR_A],sum_SHM_PCR_B$x[id_PCR_B],sum_SHM_PCR_C$x[id_PCR_C],
                                 sum_SHM_PCR_D$x[id_PCR_D],sum_SHM_PCR_E$x[id_PCR_E],sum_SHM_PCR_F$x[id_PCR_F])
colnames(reads_clones_SHM_replicates)<-c("SHM_PCR_A","SHM_PCR_B","SHM_PCR_C","SHM_PCR_D","SHM_PCR_E","SHM_PCR_F")


reads_clones_all<-cbind(reads_clones_all,reads_clones_SHM,reads_clones_SHM_isotypes,reads_clones_SHM_replicates)


###length of CDR3
cdr3_length_freq<-aggregate(data_merge$CDR3_length, by=list(data_merge$specimen_label,data_merge$amplification_template), FUN=mean)
cdr3_length_gDNA<-cdr3_length_freq[which(cdr3_length$Group.2=="gDNA"),]
cdr3_length_cDNA<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="cDNA"),]
id_gDNA<-match(rownames(reads_clones_all),cdr3_length_gDNA$Group.1)
id_cDNA<-match(rownames(reads_clones_all),cdr3_length_cDNA$Group.1)
reads_clones_CDR3<-cbind(cdr3_length_gDNA$x[id_gDNA],cdr3_length_cDNA$x[id_cDNA])
colnames(reads_clones_CDR3)<-c("CDR3_length_gDNA","CDR3_length_cDNA")

##For isotypes
cdr3_length_freq<-aggregate(data_merge$CDR3_length, by=list(data_merge$specimen_label,data_merge$isotype), FUN=mean)
cdr3_length_unmapped<-cdr3_length_freq[which(cdr3_length_freq$Group.2==""),]
cdr3_length_IGHA<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="IGHA"),]
cdr3_length_IGHD<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="IGHD"),]
cdr3_length_IGHE<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="IGHE"),]
cdr3_length_IGHG<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="IGHG"),]
cdr3_length_IGHM<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="IGHM"),]
id_unmapped<-match(rownames(reads_clones_all),cdr3_length_unmapped$Group.1)
id_IGHA<-match(rownames(reads_clones_all),cdr3_length_IGHA$Group.1)
id_IGHD<-match(rownames(reads_clones_all),cdr3_length_IGHD$Group.1)
id_IGHE<-match(rownames(reads_clones_all),cdr3_length_IGHE$Group.1)
id_IGHG<-match(rownames(reads_clones_all),cdr3_length_IGHG$Group.1)
id_IGHM<-match(rownames(reads_clones_all),cdr3_length_IGHM$Group.1)

reads_clones_CDR3_isotypes<-cbind(cdr3_length_unmapped$x[id_unmapped],cdr3_length_IGHA$x[id_IGHA],cdr3_length_IGHD$x[id_IGHD],
                                 cdr3_length_IGHE$x[id_IGHE],cdr3_length_IGHG$x[id_IGHG],cdr3_length_IGHM$x[id_IGHM])
colnames(reads_clones_CDR3_isotypes)<-c("CDR3_length_unmapped","CDR3_length_IGHA","CDR3_length_IGHD","CDR3_length_IGHE","CDR3_length_IGHG","CDR3_length_IGHM")

##For repliactes
cdr3_length_freq<-aggregate(data_merge$CDR3_length, by=list(data_merge$specimen_label,data_merge$replicates_gDNA), FUN=mean)
cdr3_length_PCR_A<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="PCR_A"),]
cdr3_length_PCR_B<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="PCR_B"),]
cdr3_length_PCR_C<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="PCR_C"),]
cdr3_length_PCR_D<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="PCR_D"),]
cdr3_length_PCR_E<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="PCR_E"),]
cdr3_length_PCR_F<-cdr3_length_freq[which(cdr3_length_freq$Group.2=="PCR_F"),]
id_PCR_A<-match(rownames(reads_clones_all),cdr3_length_PCR_A$Group.1)
id_PCR_B<-match(rownames(reads_clones_all),cdr3_length_PCR_B$Group.1)
id_PCR_C<-match(rownames(reads_clones_all),cdr3_length_PCR_C$Group.1)
id_PCR_D<-match(rownames(reads_clones_all),cdr3_length_PCR_D$Group.1)
id_PCR_E<-match(rownames(reads_clones_all),cdr3_length_PCR_E$Group.1)
id_PCR_F<-match(rownames(reads_clones_all),cdr3_length_PCR_F$Group.1)
reads_clones_CDR3_replicates<-cbind(cdr3_length_PCR_A$x[id_PCR_A],cdr3_length_PCR_B$x[id_PCR_B],cdr3_length_PCR_C$x[id_PCR_C],
                                   cdr3_length_PCR_D$x[id_PCR_D],cdr3_length_PCR_E$x[id_PCR_E],cdr3_length_PCR_F$x[id_PCR_F])
colnames(reads_clones_CDR3_replicates)<-c("CDR3_length_PCR_A","CDR3_length_PCR_B","CDR3_length_PCR_C","CDR3_length_PCR_D","CDR3_length_PCR_E","CDR3_length_PCR_F")


reads_clones_all<-cbind(reads_clones_all,reads_clones_CDR3,reads_clones_CDR3_isotypes,reads_clones_CDR3_replicates)


###To obtain the overlapping samples between clinical annotation and data
id.sample <- match(rownames(reads_clones_all),clin_annot$specimen_id)
reads_clones_annot <- cbind(clin_annot[id.sample,], reads_clones_all)
write.csv(reads_clones_annot, "/Users/Pinedasans/VDJ/Data/total_reads_clones.csv", row.names = F)

##Read the annotation file with the new subject_id
reads_clones_annot<-read.csv("/Users/Pinedasans/VDJ/Data/total_reads_clones.csv")

save(data_merge,reads_clones_annot,file="/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")

### Add naive and memory B-cells clones from the IGHM isotype
## SHM<=4 is navie SHM>=4 is memory
data_merge$IGHM_naive_memory<-ifelse(data_merge$isotype=="IGHM" & data_merge$SHM<=4,"naive",
                                     ifelse(data_merge$isotype=="IGHM" & data_merge$SHM>4,"memory",NA))

clones_naive_memory<- unique(data_merge[,c("specimen_label","V_J_lenghCDR3_CloneId","IGHM_naive_memory")])
clones_naive_memory_matrix<-data.matrix(table(clones_naive_memory$specimen_label,clones_naive_memory$IGHM_naive_memory))
colnames(clones_naive_memory_matrix)<-c("clones_memory","clones_naive")

id<-match(reads_clones_annot$specimen_id,rownames(clones_naive_memory_matrix))

reads_clones_annot$clones_naive<-clones_naive_memory_matrix[id,1]
reads_clones_annot$clones_memory<-clones_naive_memory_matrix[id,2]

sum_SHM_freq<-aggregate(data_merge$SHM_freq, by=list(data_merge$specimen_label,data_merge$IGHM_naive_memory), FUN=sum)
sum_SHM_naive<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="naive"),]
sum_SHM_memory<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="memory"),]

id_naive<-match(reads_clones_annot$specimen_id,sum_SHM_naive$Group.1)
id_memory<-match(reads_clones_annot$specimen_id,sum_SHM_memory$Group.1)

reads_clones_SHM_isotypes<-cbind(sum_SHM_naive$x[id_naive],sum_SHM_memory$x[id_memory])
colnames(reads_clones_SHM_isotypes)<-c("SHM_naive","SHM_memory")

reads_clones_annot<-cbind(reads_clones_annot,reads_clones_SHM_isotypes)

write.csv(reads_clones_annot,file="/Users/Pinedasans/VDJ/Data/total_reads_clones.csv")
save(data_merge,reads_clones_annot,file="/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")

###Put the individuals in the data_merge
load("/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")
reads_clones_annot<-read.csv("/Users/Pinedasans/VDJ/Data/total_reads_clones.csv")

xx<-match(data_merge$sample_id,reads_clones_annot$Sample_id)
data_merge$individual_id<-reads_clones_annot[xx,"Individual.id"]
save(data_merge,reads_clones_annot,file="/Users/Pinedasans/VDJ/Data/VDJ_clonesAllmerged.Rdata")


