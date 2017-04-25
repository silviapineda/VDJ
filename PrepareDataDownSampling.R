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
save(data_clonesInference, file="ClonesInfered_downsampled_amp.RData")

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

save(data_clonesInference_ind,reads_clones_down_annot,file="/Users/Pinedasans/Data/VDJ/VDJ_downSampled_amp.Rdata")

#################################################
#### ReadData DownSampling by Individual  ######
#################################################

#Read data from the output file from CloneInference.py
##Read all the files and save into and Rdata all together
files <- list.files("/Users/Pinedasans/Data/VDJ/ClonesInferedMinimum/")
data_clonesInference_ind <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.csv(paste("ClonesInferedMinimum/", i, sep = ""))
  t$dataset<-rep(i,nrow(t))
  data_clonesInference_ind <- rbind(data_clonesInference_ind, t)
}
save(data_clonesInference_ind, file="ClonesInfered_downsampled_individual.RData")

##Read counts and clones per sample and data point
read_count <- table(data_clonesInference_ind$specimen_label,data_clonesInference_ind$dataset)
read_count<-apply(read_count,1,mean)

read_count_gDNA<-NULL
read_count_cDNA<-NULL
for (i in 1:5){
  read_count_gDNA<- rbind(read_count_gDNA,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_ind$amplification_template[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  read_count_cDNA<- rbind(read_count_cDNA,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_ind$amplification_template[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
}
read_count_gDNA<-apply(read_count_gDNA,2,mean)
read_count_cDNA<-apply(read_count_cDNA,2,mean)
read_count_amplification<-cbind(read_count_gDNA,read_count_cDNA)

data_clonesInference_ind$isotype <- substr(data_clonesInference_ind$isosubtype, 1, 4)
read_count_UNMAPPED<-NULL
read_count_IGHA<-NULL
read_count_IGHD<-NULL
read_count_IGHE<-NULL
read_count_IGHG<-NULL
read_count_IGHM<-NULL
for (i in 1:5){
  read_count_UNMAPPED<- rbind(read_count_UNMAPPED,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                        data_clonesInference_ind$isotype[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,1])
  read_count_IGHA<- rbind(read_count_IGHA,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_ind$isotype[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,2])
  read_count_IGHD<- rbind(read_count_IGHD,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_ind$isotype[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,3])
  read_count_IGHE<- rbind(read_count_IGHE,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_ind$isotype[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,4])
  read_count_IGHG<- rbind(read_count_IGHG,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_ind$isotype[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,5])
  read_count_IGHM<- rbind(read_count_IGHM,table(data_clonesInference_ind$specimen_label[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))],
                                                data_clonesInference_ind$isotype[which(data_clonesInference_ind$dataset==paste("ClonesInfered_downsampled_",i,".csv",sep=""))])[,6])
  
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
data_clonesInference_ind$V_J_lenghCDR3_Clone = paste(data_clonesInference_ind$V_J_lenghCDR3,data_clonesInference_ind$numberClone,sep="_")
read_count_Clones<- unique(data_clonesInference_ind[,c("specimen_label","V_J_lenghCDR3_Clone","amplification_template","dataset")])
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
data_clonesInference_ind$V_J_lenghCDR3_Clone_igh = paste(data_clonesInference_ind$V_J_lenghCDR3,data_clonesInference_ind$igh_clone_id,sep="_")
read_count_Clones_igh<- unique(data_clonesInference_ind[,c("specimen_label","V_J_lenghCDR3_Clone_igh","amplification_template","dataset")])
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
reads_clones_down_annot_ind<- cbind(clin_annot[id.sample,], reads_clones_down)
write.csv(reads_clones_down_annot_ind, "total_reads_clones_down_ind.csv", row.names = F)

save(data_clonesInference_ind,reads_clones_down_annot_ind,file="/Users/Pinedasans/Data/VDJ/VDJ_downSampled_ind.Rdata")

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

