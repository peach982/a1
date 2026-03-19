#/Lustre01/tangqianzi/anaconda3forSeurat4/envs/R4.0.2/bin/R
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(stringr)
library(hdf5r)


args <- commandArgs()
mysample_name<- args[6]
my_input <- args[7]    #第一步中PMSfilter_f3rds文件
dimplot_path <- args[8]
spatial_dimplot_path <- args[9]
my_dim <- as.integer(args[10])  
my.R <- as.numeric(args[11])  
remove_csv <- args[12]
remove_summary_csv <- args[13]
rds_file <- args[14]


PMSfilter_f3<-readRDS(my_input)
PMSfilter_f3  <- SCTransform(PMSfilter_f3 , assay = "Spatial", return.only.var.genes = FALSE)
DefaultAssay(PMSfilter_f3) <- "SCT"
PMSfilter_f3 <- RunPCA(PMSfilter_f3)
PMSfilter_f3 <- FindNeighbors(PMSfilter_f3, dims = 1:my_dim)
PMSfilter_f3 <- FindClusters(PMSfilter_f3, resolution = my.R)
PMSfilter_f3 <- RunUMAP(PMSfilter_f3, dims = 1:my_dim)
length(unique(PMSfilter_f3$seurat_clusters))
remove_data <- read.csv(remove_csv)
rownames(remove_data) <- remove_data[,1]
PMSfilter_f3$new <- NA
PMSfilter_f3$new <- remove_data[colnames(PMSfilter_f3),2]
PMSfilter_f3 <- subset(PMSfilter_f3,subset = new %in% c('X','Y','Z','P'))
cluster_summary <- as.data.frame(summary(as.factor(remove_data[,2])))
colnames(cluster_summary) <- mysample_name
write.csv(cluster_summary, file = remove_summary_csv,quote = F, row.names = T)
p1 <- DimPlot(PMSfilter_f3, reduction = "umap",group.by = 'new',cols = c('#f58220','#33a3dc' ,'#d71345'))
p2 <- SpatialDimPlot(PMSfilter_f3,pt.size.factor = 1.55,group.by = 'new',cols =  c('#f58220','#33a3dc','#d71345'))
ggsave(filename = dimplot_path, plot=p1, width = 5, height =5)
ggsave(filename = spatial_dimplot_path, plot=p2, width = 5, height =5)
saveRDS(PMSfilter_f3, file = rds_file)
# saveRDS(PMSfilter_f3, file = paste0("/Lustre03/data/tangqianzi/Gut_timepoint_spatial/",mysample_name,"-filter+SCT+PCA.rds"))

