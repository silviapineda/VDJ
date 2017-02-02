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
### Date: November, 2016
############################################################################################
library(gplots)
require(graphics)
library(RColorBrewer)
library("bcRep")
library(ggplot2)
library(devtools)
library(easyGgplot2)


setwd("/Users/Pinedasans/Documents/VDJ/")
load("/Users/Pinedasans/Data/VDJ/VDJ.Rdata")


##Separate Long and AR
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_AR<-reads_clones_annot[which(reads_clones_annot$clin=="AR" | reads_clones_annot$clin=="pre-AR"),]

clinLong<-factor(reads_clones_annot_Long$clin, levels=c("NP", "PNR", "PR"))
clinAR<-factor(reads_clones_annot_AR$clin, levels=c("pre-AR","AR"))

####Separate by gDNA and cDNA and delete the reads < 100
reads_clones_annot_Long_cDNA <- reads_clones_annot_Long[which(reads_clones_annot_Long$cDNA_reads>=100),] ## 55 samples (Discard samples with less than 100 reads)
clinLong_cDNA <- clinLong[which(reads_clones_annot_Long$cDNA_reads>=100)] #(Discard samples with less than 100 reads)

reads_clones_annot_AR_cDNA <- reads_clones_annot_AR[which(reads_clones_annot_AR$cDNA_reads>=100),] ## 20 samples (Discard samples with less than 100 reads)
clinAR_cDNA <- clinAR[which(reads_clones_annot_AR$cDNA_reads>=100)] #(Discard samples with less than 100 reads)

reads_clones_annot_Long_gDNA <- reads_clones_annot_Long[which(reads_clones_annot_Long$gDNA_reads>=100),] ## 69 samples (Discard samples with less than 100 reads)
clinLong_gDNA <- clinLong[which(reads_clones_annot_Long$gDNA_reads>=100)] #(Discard samples with less than 100 reads)

reads_clones_annot_AR_gDNA <- reads_clones_annot_AR[which(reads_clones_annot_AR$gDNA_reads>=100),] ## 10 samples (Discard samples with less than 100 reads)
clinAR_gDNA <- clinAR[which(reads_clones_annot_AR$gDNA_reads>=100)] #(Discard samples with less than 100 reads)


###cDNA
####Plot the longitudinal data
p.cDNA.long<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=cDNA_reads, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() +
  geom_point() +
  ylim(0, 280000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("cDNA Reads")

p.cDNA.long.M154<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=M154_reads, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() +
  geom_point() +
  ylim(0, 280000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("cDNA M154 Reads")

p.cDNA.long.M155<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=M155_reads, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() +
  geom_point() +
  ylim(0, 280000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("cDNA M155 Reads")

plots <- list() 
plots[[1]] <- p.cDNA.long
plots[[2]] <- p.cDNA.long.M154
plots[[3]] <- p.cDNA.long.M155

layout <- matrix(c(1, 2, 1, 3), nrow = 2, byrow = TRUE)
ggplot2.multiplot(plotlist = plots, layout = layout)

###gDNA
p.gDNA.long<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=gDNA_reads, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("gDNA Reads")

ggplot2.multiplot(p.gDNA.long,p.cDNA.long,p.cDNA.long.M154,p.cDNA.long.M155)

###cDNA
####Plot the AR data
p.cDNA.AR<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=cDNA_reads, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  ylim(0, 280000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("cDNA Reads")

p.cDNA.AR.M154<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=M154_reads, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  ylim(0, 280000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("M154 Reads")

p.cDNA.AR.M155<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=M155_reads, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  ylim(0, 280000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("M155 Reads")

###gDNAa
p.gDNA.AR<-ggplot(data=reads_clones_annot_AR_gDNA, aes(x=time, y=gDNA_reads, group=subject_id, shape=clinAR_gDNA, color=clinAR_gDNA)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("gDNA Reads ")

# plots <- list() 
# plots[[1]] <- p.cDNA
# plots[[2]] <- p2.gDNA
# plots[[3]] <- p.gDNA
# 
# layout <- matrix(c(1, 2, 1, 3), nrow = 2, byrow = TRUE)

ggplot2.multiplot(p.gDNA.long,p.cDNA.long,p.gDNA.AR,p.cDNA.AR)


####################
####Plot Clones####

###Clones by individuals using gDNA
p1<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=clones_igh_gDNA, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("ClonesLong (per individual) gDNA")

###Clones by individuals using cDNA
p2<-ggplot(data=reads_clones_annot_Long_cDNA, aes(x=time, y=clones_igh_cDNA, group=subject_id, shape=clinLong_cDNA, color=clinLong_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("ClonesLong (per Individual cDNA")

##Clones by individuals AR using cDNA
p3<-ggplot(data=reads_clones_annot_AR_gDNA, aes(x=time, y=clones_igh_gDNA, group=subject_id, shape=clinAR_gDNA, color=clinAR_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("ClonesAR (per Individual) gDNA")

##Clones by individuals AR using cDNA
p4<-ggplot(data=reads_clones_annot_AR_cDNA, aes(x=time, y=clones_igh_cDNA, group=subject_id, shape=clinAR_cDNA, color=clinAR_cDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("goldenrod","firebrick3")) +
  ggtitle("ClonesAR (per Individual) cDNA")

ggplot2.multiplot(p1,p2,p3,p4)


##Clones by samples using gDNA
p1.sample<-ggplot(data=reads_clones_annot_Long_gDNA, aes(x=time, y=total_clones_gDNA, group=subject_id, shape=clinLong_gDNA, color=clinLong_gDNA)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("ClonesLong (per sample) gDNA")




################################
##### Demographics for table ###
###############################
table(reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

##Recipient Age
summary(reads_clones_annot$Recipient.Age.when.had.Tx[!duplicated(reads_clones_annot$Sample_id)])
#NP
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
summary(reads_clones_annot$Recipient.Age.when.had.Tx[which(reads_clones_annot$clin=="AR")])

fit = lm(reads_clones_annot$Recipient.Age.when.had.Tx[!duplicated(reads_clones_annot$Sample_id)] ~ reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])
anova(fit)



##Donor Age
summary(reads_clones_annot$Donor.Age[!duplicated(reads_clones_annot$Sample_id)])
#NP
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
summary(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin=="AR")])

fit = lm(reads_clones_annot$Donor.Age[!duplicated(reads_clones_annot$Sample_id)] ~ reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])
anova(fit)

fit = lm(reads_clones_annot$Donor.Age[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)] ~
      reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])
anova(fit)

##immunosuppression
table(reads_clones_annot$immunosuppression[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$immunosuppression[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$immunosuppression[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
      reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])

##Recipient.Gender
table(reads_clones_annot$Recipient.Gender[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$Recipient.Gender[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$Recipient.Gender[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])


##Donor.Gender
table(reads_clones_annot$Donor.Gender[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$Donor.Gender[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$Donor.Gender[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])


##recipient.Race
table(reads_clones_annot$recipient.Race[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$recipient.Race[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$recipient.Race[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])

##Relatedness
table(reads_clones_annot$Donor.Source[!duplicated(reads_clones_annot$Sample_id)])
#NP
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="NP")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="NP")]])
#PNR
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="PNR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PNR")]])
#PR
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="PR")][!duplicated(reads_clones_annot$Sample_id)[which(reads_clones_annot$clin=="PR")]])
#AR
table(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin=="AR")])

chisq.test(reads_clones_annot$Donor.Source[!duplicated(reads_clones_annot$Sample_id)], reads_clones_annot$clin[!duplicated(reads_clones_annot$Sample_id)])

chisq.test(reads_clones_annot$Donor.Source[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)], 
           reads_clones_annot$clin[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR")][!duplicated(reads_clones_annot$Sample_id)])



###Considering the number of clones gDNA
par(mfrow = c(2, 2))  
boxplot(reads_clones_annot_Long_gDNA$clones_igh[which(reads_clones_annot_Long_gDNA$time==0)] ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)], 
        main="Time 0", ylab = "Number of clones", col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_gDNA$clones_igh[which(reads_clones_annot_Long_gDNA$time==0)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==0)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clones_igh[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)], 
            main="Time 6",ylab = "Number of clones", col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_gDNA$clones_igh[which(reads_clones_annot_Long_gDNA$time==6)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==6)])
anova(fit)

boxplot((reads_clones_annot_Long_gDNA$clones_igh[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)], ylab = "Number of clones", 
            main="Time 24",col = c("chartreuse4", "dodgerblue3","darkorange2"))
fit = lm((reads_clones_annot_Long_gDNA$clones_igh[which(reads_clones_annot_Long_gDNA$time==24)]) ~ clinLong_gDNA[which(reads_clones_annot_Long_gDNA$time==24)])
anova(fit)

boxplot(reads_clones_annot_AR_gDNA$clones_igh ~ clinAR_gDNA, main = "AR" , ylab = "Number of clones",col = c("goldenrod","firebrick3"))
fit = lm((reads_clones_annot_AR_gDNA$clones_igh~ clinAR_gDNA))
anova(fit)

###Considering the number of clones cDNA
par(mfrow = c(1, 2))  
boxplot(reads_clones_annot_AR_gDNA$clones_igh ~ clinAR_gDNA, main = "AR gDNA" , ylab = "Number of clones",col = c("goldenrod","firebrick3"))
boxplot(reads_clones_annot_AR_cDNA$clones_igh ~ clinAR_cDNA, main = "AR cDNA" , ylab = "Number of clones",col = c("goldenrod","firebrick3"))
fit = lm((reads_clones_annot_AR_cDNA$clones_igh~ clinAR_cDNA))
anova(fit)




