library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(stringr)
library(hdf5r)

args <- commandArgs()
n <- length(args)

mysample_name<- args[6]    
my_input <- args[7]
my_dim <- as.integer(args[8])  
my.R <- as.numeric(args[9])  
summary_csv_path <- args[10]
dimplot_path <- args[11]
spatial_dimplot_path <- args[12]
cluster_csv_path <- args[13]
rds_path <- args[14]
my_bin <- args[15:n]

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

#p1 <- DimPlot(PMSfilter_f3, reduction = "umap", group.by = 'seurat_clusters')
#p2 <- SpatialDimPlot(PMSfilter_f3,pt.size.factor = 1.55,group.by = 'seurat_clusters')
#ggsave(filename = dimplot_path, plot=p1, width = 5, height =5)
#ggsave(filename = spatial_dimplot_path, plot=p2, width = 5, height =5)




XYZP.cluster <- as.vector(PMSfilter_f3$seurat_clusters)
n <- length(my_bin) / 2
for (i in 1:n) {
	XYZP.cluster[XYZP.cluster== as.integer(my_bin[(i*2)-1])] <- my_bin[i*2]
}

#XYZP.cluster[XYZP.cluster== 3] <- "UN"

cluster_summary <- as.data.frame(summary(as.factor(XYZP.cluster)))
colnames(cluster_summary) <- mysample_name
write.csv(cluster_summary, file = summary_csv_path,quote = F, row.names = T)

#write.csv(cluster_summary, file = paste0("/Lustre03/data/tangqianzi/Gut_timepoint_spatial/",mysample_name,"_XYZP_summary.csv"),quote = F, row.names = T)

PMSfilter_f3$new <- as.factor(XYZP.cluster)
cluster.info <- data.frame(PMSfilter_f3$new)

p1 <- DimPlot(PMSfilter_f3, reduction = "umap",group.by = 'new',cols = c('purple','#f58220','#33a3dc' ,'#d71345'))
p2 <- SpatialDimPlot(PMSfilter_f3,pt.size.factor = 1.55,group.by = 'new',cols =  c('purple','#f58220','#33a3dc','#d71345'))

ggsave(filename = dimplot_path, plot=p1, width = 5, height =5)
ggsave(filename = spatial_dimplot_path, plot=p2, width = 5, height =5)

#ggsave(filename = paste0("/Lustre03/data/tangqianzi/Gut_timepoint_spatial/",mysample_name,"_final_DimPlot.pdf"), plot=p1, width = 5, height =5)
#ggsave(filename = paste0("/Lustre03/data/tangqianzi/Gut_timepoint_spatial/",mysample_name,"_final_SpatialDimPlot.pdf"), plot=p2, width = 5, height =5)


colnames(cluster.info) <- paste0(mysample_name,"_XYZ")
colnames(PMSfilter_f3@meta.data)[grep("new",colnames(PMSfilter_f3@meta.data))] <- paste0(mysample_name,"_XYZ")
write.csv(cluster.info, file = cluster_csv_path,quote = F, row.names = T)

#write.csv(cluster.info, file = paste0("/Lustre03/data/tangqianzi/Gut_timepoint_spatial/",mysample_name,"_XYZP_cluster.csv"),quote = F, row.names = T)
#p3 <- p1 + p2
#plot.name <- paste0("plot/",mysample_name,"-cluster-R",my.R,"-dims10-UMAP_and_SpatialPlot.pdf")
#ggsave(filename = plot.name, plot=p3, width = 12, height =5)
saveRDS(PMSfilter_f3, file = rds_path)


#saveRDS(PMSfilter_f3, file = paste0("/Lustre03/data/tangqianzi/Gut_timepoint_spatial/",mysample_name,"-filter+SCT+PCA.rds"))
