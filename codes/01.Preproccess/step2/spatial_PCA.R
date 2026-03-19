library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(stringr)
library(hdf5r)

args <- commandArgs()

mysample_name<- args[6]    
my_input <- args[7]
my_output <- args[8]

#rawinfo <- PMSfilter@meta.data[,c(1:4,18)]
#rawinfo$samInfo <- paste0("raw_",rawinfo$samInfo)
#cleaninfo <- PMSfilter_f3@meta.data[,c(1:4,18)]
#cleaninfo$samInfo <- paste0("clean_",cleaninfo$samInfo)
#write.csv(rawinfo, file = paste0("sampleInfo/",mysample_name,"-spotsInfo-raw.csv"))
#write.csv(cleaninfo, file = paste0("sampleInfo/",mysample_name,"-spotsInfo-filter.csv"))

#PMSfilter
#PMSfilter_f1
#PMSfilter_f2
#PMSfilter_f3

#saveRDS(PMSfilter, file = paste0("rawdata/",mysample_name,".rds"))

PMSfilter_f3<-readRDS(my_input)
PMSfilter_f3  <- SCTransform(PMSfilter_f3 , assay = "Spatial", return.only.var.genes = FALSE)
DefaultAssay(PMSfilter_f3) <- "SCT"

PMSfilter_f3 <- RunPCA(PMSfilter_f3)
p1 <- ElbowPlot(PMSfilter_f3)
ggsave(filename = my_output, plot=p1, width = 5, height =5)
