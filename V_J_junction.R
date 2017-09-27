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
### Date: June, 2017
############################################################################################
library(circlize)
library("RColorBrewer")
library(gtools)

setwd("/Users/Pinedasans/VDJ/Results/")
load("/Users/Pinedasans/VDJ/Data/VDJ_order.Rdata")

############
### gDNA ###
############
data_qc_gDNA<-data_qc[which(data_qc$amplification_template=="gDNA"),]
### longitudinal
data_qc_gDNA_long<-data_qc_gDNA[which(data_qc_gDNA$clin!="AR" & data_qc_gDNA$clin!="pre-AR"),]
reads_clones_annot_Long<-reads_clones_annot[which(reads_clones_annot$clin!="AR" & reads_clones_annot$clin!="pre-AR"),]
reads_clones_annot_Long_qc<-reads_clones_annot_Long[which(reads_clones_annot_Long$reads_gDNA>100),]

id<-match(data_qc_gDNA_long$specimen_label,reads_clones_annot_Long_qc$specimen_id)
data_qc_gDNA_long_qc<-data_qc_gDNA_long[which(is.na(id)==F),]
data_qc_gDNA_long_qc$specimen_label<-factor(data_qc_gDNA_long_qc$specimen_label)


###Circular plot for V and J junctionfor(i in 1:length(specimen)){
specimen<-unique(data_qc_gDNA_long_qc$specimen_label)
for(i in 1:length(specimen)){
  print(i)
  data_qc_gDNA_long_qc_specimen<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label==specimen[i]),]
  v_j_gene_mat<-unclass(table(data_qc_gDNA_long_qc_specimen$v_gene,data_qc_gDNA_long_qc_specimen$j_gene))
  circos.clear()
  tiff(paste("v_j_junction_gDNA_long_",specimen[i],".tiff",sep=""),res=300,h=1800,w=1800)
  circos.par(gap.degree = c(rep(2, nrow(v_j_gene_mat)-1), 10,
                          rep(2, ncol(v_j_gene_mat)-1), 10))

  grid.col<-c(rev(colorRampPalette(brewer.pal(9,"Set1"))(nrow(v_j_gene_mat))),
            colorRampPalette(brewer.pal(6,"Blues"))(ncol(v_j_gene_mat)))

  names(grid.col)<-c(as.character(rownames(v_j_gene_mat)),as.character(colnames(v_j_gene_mat)))

  chordDiagram(v_j_gene_mat, grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))
  # we go back to the first track and customize sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
  }, bg.border = NA) 

  dev.off()
} 



###Using data frames
specimen<-unique(data_qc_gDNA_long_qc$specimen_label)
for(i in 1:length(specimen)){
  print(i)
  data_qc_gDNA_long_qc_specimen<-data_qc_gDNA_long_qc[which(data_qc_gDNA_long_qc$specimen_label==specimen[i]),]
  v_j_gene_df<-data.frame(table(data_qc_gDNA_long_qc_specimen$v_gene,data_qc_gDNA_long_qc_specimen$j_gene))
  v_j_gene_df.nonceros<-v_j_gene_df[which(v_j_gene_df$Freq!=0),]
  v_j_gene_df.nonceros_oder<-v_j_gene_df.nonceros[order(v_j_gene_df.nonceros[,1]), ]
  circos.clear()
  tiff(paste("v_j_junction_gDNA_long_",specimen[i],".tiff",sep=""),res=300,h=1800,w=1800)
  circos.par(gap.degree = c(rep(2, length(unique(v_j_gene_df.nonceros_oder[[1]]))-1), 10,
                          rep(2, length(unique(v_j_gene_df.nonceros_oder[[2]]))-1), 10))

  grid.col<-c(colorRampPalette(brewer.pal(9,"Set1"))(length(unique(v_j_gene_df.nonceros_oder[[1]]))),
              colorRampPalette(brewer.pal(6,"Blues"))(length(unique(v_j_gene_df.nonceros_oder[[2]]))))
  names(grid.col)<-c(as.character(unique(v_j_gene_df.nonceros_oder[[1]])),as.character(unique(v_j_gene_df.nonceros_oder[[2]])))

  chordDiagram(v_j_gene_df.nonceros_oder, grid.col = grid.col, annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))
  # we go back to the first track and customize sector labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5),cex = 0.6)
  }, bg.border = NA) 
  dev.off()
} 



