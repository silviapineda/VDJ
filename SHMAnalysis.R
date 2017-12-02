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
### DESCRIP: SHM analysis
###         
###
### Author: Silvia Pineda
### Date: November, 2017
############################################################################################

setwd("/Users/Pinedasans/VDJ/SummaryResults/Vgene/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")

############
### gDNA ###
############
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
