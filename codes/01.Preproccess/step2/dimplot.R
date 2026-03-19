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

my_dim <- as.integer(args[8])  
my.R <- as.numeric(args[9])  
dimplot_path <- args[10]
spatial_dimplot_path <- args[11]


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
#p1 <- ElbowPlot(PMSfilter_f3)
#ggsave(filename = my_output, plot=p1, width = 5, height =5)


##p1 <- DimHeatmap(PMSfilter_f3, dims = c(1,2,3,4,5,6), cells = 2000, balanced = TRUE,fast = F)
##ggsave(filename = paste0("plot/",mysample_name,"_DimHeatmap.pdf"), plot=p1, width = 18, height =8)
PMSfilter_f3 <- FindNeighbors(PMSfilter_f3, dims = 1:my_dim)

#}

##PMSfilter_f3 <- readRDS("rawdata/sc3R-filter+SCT+PCA.rds")

PMSfilter_f3 <- FindClusters(PMSfilter_f3, resolution = my.R)
PMSfilter_f3 <- RunUMAP(PMSfilter_f3, dims = 1:my_dim)

length(unique(PMSfilter_f3$seurat_clusters))

p1 <- DimPlot(PMSfilter_f3, reduction = "umap", group.by = 'seurat_clusters')
p2 <- SpatialDimPlot(PMSfilter_f3,pt.size.factor = 1.55,group.by = 'seurat_clusters')
ggsave(filename = dimplot_path, plot=p1, width = 5, height =5)
ggsave(filename = spatial_dimplot_path, plot=p2, width = 5, height =5)