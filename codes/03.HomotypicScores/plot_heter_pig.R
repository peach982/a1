
library(Seurat)
##library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(imager)
library(EBImage)
########pig
sample<-c("1D","2J_1","2J_2","3I_1","3I_2","4CA_1","4CA_2","4CA_3","5CM","6PC_1","6PC_2","6PC_3","7DC_1")
colors <- c(
  "Enterocyte" = "#c99fc7",
  "Goblet" = "#ece863",
  "Stem" = "#a2d19c",
  "EEC" = "#4b6aa8",  # 使用最后出现的EEC颜色
  "Myeloid" = "#e8905e",
  "T_NKcell" = "#53738c",
  "Bcell" = "#7ec9c4",
  "Plasma" = "#408444",
  "TA" = "#b7deea",
  "Paneth" = "#6a9955",
  "Neutrophil" = "#519aba",
  "SecretoryProgenitors" = "#aa4136",
  "InnateLymphoid" = "#d25774",
  "AdaptiveLymphoid" = "#e6e2a3"
)

for (m in sample){
file_path <- paste0('/data/liuqinglin/ST/spatial/data_spatial/rds/', m, "-filter+SCT+PCA.rds")
brain <- readRDS(file_path)
cluster_meta<-read.table(file=paste0('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_metadata/test_',m,'_metadata.subset.xls'),sep='\t',header=T,row.names=1)
mykeepcells<-rownames(cluster_meta)
library(stringr)
mykeepcells2<-as.vector(str_split_fixed(mykeepcells,'_',2)[,1])
brain2 <-subset(brain,cells=mykeepcells2)
check1<-brain2@images$PMSfilter@image
dim(check1)<-c(dim(check1)[1],dim(check1)[2],1,3)
check<-grayscale(check1)
dim(check)<-c(dim(check1)[1],dim(check1)[2])

cells = rgbImage(check, check, check)
cells2 = aperm(cells, c(2,1,3))
cells2 = aperm(cells2, c(2,1,3))
brain2@images$PMSfilter@image<-cells2

keepcelltype<-rownames(cluster_meta[cluster_meta$seurat_clusters%in%celltype,])
library(stringr)
keepcelltype2<-as.vector(str_split_fixed(keepcelltype,'_',2)[,1])
brain_celltype <-subset(brain2,cells=keepcelltype2)

mydegreenames<-rownames(cluster_meta)
library(stringr)
mydegreenames2<-as.vector(str_split_fixed(mydegreenames,'_',2)[,1])
rownames(cluster_meta)<-mydegreenames2

brain_celltype@meta.data$celltype<-as.vector(cluster_meta[rownames(brain_celltype@meta.data),6])
###Enterocytes #FA8072#  Goblet #87CEFA#  EECs #D28EFF# TA #32CD32# Progenior  #FFFACD#  Stem.Cell #FFC0CB#
library(RColorBrewer)


Idents(brain_celltype) <- brain_celltype@meta.data$celltype

p1<-SpatialPlot(brain_celltype, group.by = "celltype", crop=FALSE, pt.size.factor = 1,cols=colors) + theme(legend.position = "right") 
ggsave(filename = paste0('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_metadata/heter_scores/',m,'_','2.pdf'), plot=p1, width = 6, height =6)
}
