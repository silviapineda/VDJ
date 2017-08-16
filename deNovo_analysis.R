rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies. deNovo analysis
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

##############################
### gDNA longitudinal data ###
##############################

data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
#data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

### longitudinal
data_qc_gDNA_long<-data_qc_gDNA[which(data_qc_gDNA$clin!="AR" & data_qc_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_gDNA>100),]

id<-match(data_qc_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_gDNA_long_qc<-data_qc_gDNA_long[which(is.na(id)==F),]
data_qc_gDNA_long_qc$specimen_label<-factor(data_qc_gDNA_long_qc$specimen_label)

###Delete the samples that do not have three time points

###7_s63
data_qc_gDNA_long_qc2<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label!="7_S63"),]
##38021
data_qc_gDNA_long_qc3<-data_qc_gDNA_long_qc2[which(data_qc_gDNA_long_qc2$sample_id!="38021" & data_qc_gDNA_long_qc2$sample_id!="12008"
                                                   & data_qc_gDNA_long_qc2$sample_id!="38020" & data_qc_gDNA_long_qc2$sample_id!="41006"
                                                   & data_qc_gDNA_long_qc2$sample_id!="89002" & data_qc_gDNA_long_qc2$sample_id!="12001"
                                                   & data_qc_gDNA_long_qc2$sample_id!="89015" & data_qc_gDNA_long_qc2$sample_id!="301003"
                                                   & data_qc_gDNA_long_qc2$sample_id!="303001"),]

data_qc_gDNA_long_clean<-data_qc_gDNA_long_qc3

##Obtain deNovo number of clones
sample_id<-unique(data_qc_gDNA_long_clean$sample_id)
deNovo<-NULL
persistance<-NULL
deNovo6<-NULL
persistance6<-NULL
deNovo24<-NULL
persistance24<-NULL
for (i in 1:length(sample_id)){
  print(i)
  clone_type_matrix<-as.data.frame.matrix(table(data_qc_gDNA_long_clean$V_J_lenghCDR3_Clone_igh[which(data_qc_gDNA_long_clean$sample_id==sample_id[i])],
                                      factor(data_qc_gDNA_long_clean$specimen_label[which(data_qc_gDNA_long_clean$sample_id==sample_id[i])])))
  clone_type_matrix$deNovo<-ifelse(clone_type_matrix[,1] == 0 & (clone_type_matrix[,2] !=0 | clone_type_matrix[,3] !=0),1,0)
  clone_type_matrix$persistance24<-ifelse(clone_type_matrix[,1] > 0 & clone_type_matrix[,2] > 0 & clone_type_matrix[,3] >0,1,0)
  clone_type_matrix$deNovo6<-ifelse(clone_type_matrix[,1] == 0 & clone_type_matrix[,2] !=0,1,0)
  clone_type_matrix$deNovo24<-ifelse(clone_type_matrix[,1] == 0  & clone_type_matrix[,2] !=0 & clone_type_matrix[,3] !=0,1,0)
  clone_type_matrix$persistance6<-ifelse(clone_type_matrix[,1] > 0 & clone_type_matrix[,2] > 0 ,1,0)
  clone_type_matrix$persistance<-ifelse(clone_type_matrix[,1] > 0 & (clone_type_matrix[,2] > 0 | clone_type_matrix[,3] > 0),1,0)
  deNovo[i]<-sum(clone_type_matrix$deNovo)
  persistance[i]<-sum(clone_type_matrix$persistance)
  deNovo6[i]<-sum(clone_type_matrix$deNovo6)
  persistance6[i]<-sum(clone_type_matrix$persistance6)
  deNovo24[i]<-sum(clone_type_matrix$deNovo24)
  persistance24[i]<-sum(clone_type_matrix$persistance24)
}
results<-cbind(deNovo,deNovo6,deNovo24,persistance,persistance6,persistance24)
rownames(results)<-sample_id
results<-data.frame(results)


###Read the clean annotated data
annot<-read.csv("/Users/Pinedasans/VDJ/Data/Annot_clean_long_gDNA.csv")
data_annot<-annot[match(rownames(results),annot$Sample_id),c(1,2,4,6:24)]

results$clin<-data_annot$clin
cols<-c("chartreuse4", "dodgerblue3","darkorange2")[results$clin]
tiff("deNovo_persistence_clin_gDNA.tiff",h=3000,w=2000,res=300)
par(mfrow = c(3,2))
barplot(results$deNovo,col=cols,main="deNovo",xlab = "Samples", ylab = "Number of deNovo Clones",las=2,
        names.arg=rownames(results))
legend(18, 6500, legend=levels(results$clin),col=c("chartreuse4", "dodgerblue3","darkorange2"),pch=15, cex=0.8)
barplot(results$persistance,col=cols,main="persistence", ylab = "Number of deNovo Clones",las=2,
        names.arg=rownames(results))
barplot(results$deNovo6,col=cols,main="deNovo6", ylab = "Number of deNovo Clones",las=2,
        names.arg=rownames(results))
barplot(results$persistance6,col=cols,main="persistence6", ylab = "Number of deNovo Clones",las=2,
        names.arg=rownames(results))
barplot(results$deNovo24,col=cols,main="deNovo24", ylab = "Number of deNovo Clones",las=2,
        names.arg=rownames(results))
barplot(results$persistance24,col=cols,main="persistence24",ylab = "Number of deNovo Clones",las=2,
        names.arg=rownames(results))
dev.off()

#clin
summary(glm(deNovo~clin,data=results))
anova(lm(deNovo~clin,data=results))

summary(glm(deNovo6~clin,data=results))
anova(lm(deNovo6~clin,data=results))

summary(glm(deNovo24~clin,data=results))
anova(lm(deNovo24~clin,data=results))

summary(glm(persistance~clin,data=results))
anova(lm(persistance~clin,data=results))

summary(glm(persistance6~clin,data=results))
anova(lm(persistance~clin,data=results))

summary(glm(persistance24~clin,data=results))
anova(lm(persistance~clin,data=results))


tiff("boxplot_deNovopersistance_clin_gDNA.tiff",h=3000,w=2000,res=300)
par(mfrow = c(3,2))
boxplot(deNovo~clin,data=results,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="deNovo")
boxplot(persistance~clin,data=results,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="persistance")
boxplot(deNovo6~clin,data=results,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="deNovo6")
boxplot(persistance6~clin,data=results,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="persistance6")
boxplot(deNovo24~clin,data=results,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="deNovo24")
boxplot(persistance24~clin,data=results,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="persistance24")
dev.off()


##Multivariable analysis
library(MASS)
COLOR=brewer.pal(5,"Set1")

####deNovo
fit<-lm(deNovo~clin+immunosuppression+GenderMismatch+recipient.Race+
              Recipient.Age.when.had.Tx+Donor.Source+hla_mismatch,data=data_annot)
step<-stepAIC(fit,direction="both")
step$anova
#Final Model:
#deNovo ~ GenderMismatch + Recipient.Age.when.had.Tx

summary(glm(deNovo~immunosuppression + GenderMismatch + Recipient.Age.when.had.Tx + 
              Donor.Source, data = data_annot))

tiff("boxplot_deNovo_GenderMismatch_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$GenderMismatch,col=COLOR,main="deNovo Clones by GenderMismatch")
dev.off()

tiff("boxplot_Recipient.Age.when.had.Tx_gDNA.tiff",h=2000,w=1800,res=300)
plot(deNovo~data_annot$Recipient.Age.when.had.Tx,pch = 16, cex = 1.3, col = "blue",main="deNovo Clones by Recipient Race")
abline(lm(deNovo~data_annot$Recipient.Age.when.had.Tx))
dev.off()

###deNovo6
fit<-lm(deNovo6~clin+immunosuppression+GenderMismatch+
          Recipient.Age.when.had.Tx+Donor.Source+hla_mismatch,data=data_annot)
step<-stepAIC(fit,direction="both")
step$anova
#Final Model:
#deNovo6 ~ immunosuppression + GenderMismatch + Recipient.Age.when.had.Tx + hla_mismatch

summary(glm(deNovo6~immunosuppression + GenderMismatch + Recipient.Age.when.had.Tx + hla_mismatch, data = data_annot))

tiff("boxplot_deNovo6_immunosuppression_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$immunosuppression,col=COLOR,main="deNovo6 Clones by immunosuppression")
dev.off()

tiff("boxplot_deNovo6_GenderMismatch_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$GenderMismatch,col=COLOR,main="deNovo6 Clones by GenderMismatch")
dev.off()

tiff("boxplot_deNovo6_Recipient.Age.when.had.Tx_gDNA.tiff",h=2000,w=1800,res=300)
plot(deNovo6~data_annot$Recipient.Age.when.had.Tx,pch = 16, cex = 1.3, col = "blue",main="deNovo6 Clones by Recipient Race")
abline(lm(deNovo6~data_annot$Recipient.Age.when.had.Tx))
dev.off()

tiff("boxplot_deNovo6_hla_mismatch_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$hla_mismatch,col=COLOR,main="deNovo6 Clones by hla_mismatch")
dev.off()

#deNovo24
fit<-lm(deNovo24~clin+immunosuppression+GenderMismatch+
          Recipient.Age.when.had.Tx+Donor.Source+hla_mismatch,data=data_annot)
step<-stepAIC(fit,direction="both")
step$anova
#Final Model:
#  deNovo24 ~ immunosuppression + hla_mismatch
summary(glm(deNovo24~immunosuppression + hla_mismatch, data = data_annot))

##persistance
fit<-lm(persistance~immunosuppression+GenderMismatch+
          Recipient.Age.when.had.Tx+Donor.Source+hla_mismatch,data=data_annot)
step<-stepAIC(fit,direction="both")
step$anova
#persistance ~ clin + GenderMismatch + Recipient.Age.when.had.Tx + Donor.Source
summary(glm(persistance~GenderMismatch + Recipient.Age.when.had.Tx + 
              Donor.Source, data = data_annot))
tiff("boxplot_persistance_genderMismatch_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(persistance~data_annot$GenderMismatch,col=COLOR,main="persistence Clones by GenderMismatch")
dev.off()

tiff("boxplot_persistence_Recipient.Age.when.had.Tx_gDNA.tiff",h=2000,w=1800,res=300)
plot(persistance~data_annot$Recipient.Age.when.had.Tx,pch = 16, cex = 1.3, col = "blue",main="persistence Clones by Recipient Race")
abline(lm(persistance~data_annot$Recipient.Age.when.had.Tx))
dev.off()

tiff("boxplot_persistance_Donor.Source_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(persistance~data_annot$Donor.Source,col=COLOR,main="persistence Clones by Donor.Source")
dev.off()

###persistence6
fit<-lm(persistance6~clin+immunosuppression+GenderMismatch+
          Recipient.Age.when.had.Tx+Donor.Source+hla_mismatch,data=data_annot)
step<-stepAIC(fit,direction="both")
step$anova
#Final Model:
 # persistance6 ~ immunosuppression + GenderMismatch + Recipient.Age.when.had.Tx + Donor.Source + hla_mismatch
summary(glm(persistance6~immunosuppression +  GenderMismatch + Recipient.Age.when.had.Tx + Donor.Source, data = data_annot))

tiff("boxplot_persistance6_immunosuppression_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(persistance6~data_annot$immunosuppression,col=COLOR,main="persistence6 Clones by immunosuppression")
dev.off()

tiff("boxplot_persistance6_genderMismatch_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(persistance6~data_annot$GenderMismatch,col=COLOR,main="persistence6 Clones by GenderMismatch")
dev.off()

tiff("boxplot_persistence6_Recipient.Age.when.had.Tx_gDNA.tiff",h=2000,w=1800,res=300)
plot(persistance6~data_annot$Recipient.Age.when.had.Tx,pch = 16, cex = 1.3, col = "blue",main="persistence6 Clones by Recipient Race")
abline(lm(persistance6~data_annot$Recipient.Age.when.had.Tx))
dev.off()

tiff("boxplot_persistance6_Donor.Source_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(persistance6~data_annot$Donor.Source,col=COLOR,main="persistence6 Clones by Donor.Source")
dev.off()


####Indiviudal variables

#immunosuppression
COLOR=brewer.pal(5,"Set1")
tiff("boxplot_deNovo_immuno_gDNA.tiff",h=2000,w=1800,res=300)
summary(glm(persistance24~data_annot$immunosuppression))
boxplot(persistance24~data_annot$immunosuppression,col=COLOR,main="deNovo Clones by Immunosuppression")
dev.off()

#Recipient Race
summary(glm(deNovo24~data_annot$recipient.Race))
tiff("boxplot_deNovo_recRace_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo24~data_annot$recipient.Race,col=COLOR,main="deNovo Clones by Recipient Race")
dev.off()

#Gender mismatch
summary(glm(deNovo~data_annot$GenderMismatch))
tiff("boxplot_deNovo_GenderMismatch_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$GenderMismatch,col=COLOR,main="deNovo Clones by Gender Mismatch")
dev.off()

#Age Recipient
summary(glm(deNovo~data_annot$Recipient.Age.when.had.Tx))
tiff("boxplot_deNovo_AgeRec_gDNA.tiff",h=2000,w=1800,res=300)
plot(deNovo~data_annot$Recipient.Age.when.had.Tx,pch = 16, cex = 1.3, col = "blue",main="deNovo Clones by Recipient Race")
abline(lm(deNovo~data_annot$Recipient.Age.when.had.Tx))
dev.off()

#Donor.Source
summary(glm(deNovo~data_annot$Donor.Source))
tiff("boxplot_deNovo_DonorSource_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$Donor.Source,col=COLOR,main="deNovo Clones by Donor Source")
dev.off()

#hla_mismatch
summary(glm(deNovo~data_annot$hla_mismatch))
tiff("boxplot_deNovo_hla_gDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$hla_mismatch,col=COLOR,main="deNovo Clones by HLA mismatch")
dev.off()



##############################
### cDNA longitudinal data ###
##############################

data_qc_cDNA<-data_qc[which(data_qc$amplification_template=="cDNA"),]

### longitudinal
data_qc_cDNA_long<-data_qc_cDNA[which(data_qc_cDNA$clin!="AR" & data_qc_cDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_cDNA>100),]

id<-match(data_qc_cDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_cDNA_long_qc<-data_qc_cDNA_long[which(is.na(id)==F),]
data_qc_cDNA_long_qc$specimen_label<-factor(data_qc_cDNA_long_qc$specimen_label)

###Delete the samples that do not have two time points

###
data_qc_cDNA_long_qc<-data_qc_cDNA_long_qc[which(data_qc_cDNA_long_qc$specimen_label!="7_S31"),]
data_qc_cDNA_long_qc2<-data_qc_cDNA_long_qc[which(data_qc_cDNA_long_qc$specimen_label!="7_S63"),]

##38021
data_qc_cDNA_long_qc3<-data_qc_cDNA_long_qc2[which(data_qc_cDNA_long_qc2$sample_id!="1005" & data_qc_cDNA_long_qc2$sample_id!="303010"
                                                   & data_qc_cDNA_long_qc2$sample_id!="38020"),]

data_qc_cDNA_long_clean<-data_qc_cDNA_long_qc3

##Obtain deNovo number of clones
sample_id<-unique(data_qc_cDNA_long_clean$sample_id)
deNovo<-NULL
expansion<-NULL
for (i in 1:length(sample_id)){
  print(i)
  clone_type_matrix<-as.data.frame.matrix(table(data_qc_cDNA_long_clean$V_J_lenghCDR3_Clone_igh[which(data_qc_cDNA_long_clean$sample_id==sample_id[i])],
                                                factor(data_qc_cDNA_long_clean$specimen_label[which(data_qc_cDNA_long_clean$sample_id==sample_id[i])])))
  clone_type_matrix$deNovo<-ifelse(clone_type_matrix[,1] == 0 & (clone_type_matrix[,2] !=0),1,0)
  clone_type_matrix$expansion<-ifelse(clone_type_matrix[,1] > 0 & clone_type_matrix[,2] >0,1,0)
  deNovo[i]<-sum(clone_type_matrix$deNovo)
  expansion[i]<-sum(clone_type_matrix$expansion)
}
names(deNovo)<-sample_id
names(expansion)<-sample_id


###Read the clean annotated data
annot<-read.csv("/Users/Pinedasans/VDJ/Data/Annot_clean_long_gDNA.csv")
data_annot<-annot[match(names(deNovo),annot$Sample_id),c(1,2,4,6:24)]

#clin deNovo
summary(glm(deNovo~data_annot$clin))
anova(lm(deNovo~data_annot$clin))
tiff("boxplot_deNovo_clin_cDNA.tiff",h=2000,w=1800,res=300)
boxplot(deNovo~data_annot$clin,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="deNovo Clones by Clin")
dev.off()

#clin expansion
summary(glm(expansion~data_annot$clin))
anova(lm(expansion~data_annot$clin))
tiff("boxplot_expansion_clin_cDNA.tiff",h=2000,w=1800,res=300)
boxplot(expansion~data_annot$clin,col=c("chartreuse4", "dodgerblue3","darkorange2"),main="Expansion Clones by Clin")
dev.off()
