library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(stringr)
library(hdf5r)

#list <- list.files("D:/ST_KFL_lxh_20210819/")
#list

args <- commandArgs()

my_dat <- data.frame()   

mysample_name <- args[6]    ########TODO
my_path <- args[7]
res_path <- args[8]
#{

#data.dir = paste0("D:/ST_KFL_lxh_20210819/", mysample_name)
data.dir = my_path
setwd(data.dir)

# dir.create("rawdata")
# dir.create("plot")
# dir.create("sampleInfo")
#
getwd()
PMSfilter <- Load10X_Spatial("outs", assay="Spatial", slice="PMSfilter", filter.matrix=TRUE)

my_dat[1,1] <- 1      #TODO
my_dat[1,2] <- dim(PMSfilter)[1]    #TODO   
my_dat[1,3] <- dim(PMSfilter)[2]   #TODO

#rawcounts <- as.data.frame(as.matrix(PMSfilter@assays$Spatial@counts))
#counts.file <- paste0("rawdata/",mysample_name, "rawcounts.csv")
#write.csv(rawcounts, file=counts.file)


#summarize spot number, total gene number, gene number per spot before filter

samInfo <- rep(mysample_name, nrow(PMSfilter@meta.data))
PMSfilter <- AddMetaData(object = PMSfilter, metadata = samInfo, col.name = "samInfo")

PMSfilter[["ND6"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ND6")
PMSfilter[["ND5"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ND5")
PMSfilter[["ND4L"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ND4L")
PMSfilter[["ND4"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ND4")
PMSfilter[["ND3"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ND3")
PMSfilter[["ND2"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ND2")
PMSfilter[["ND1"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ND1")
PMSfilter[["CYTB"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "CYTB")
PMSfilter[["COX3"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "COX3")
PMSfilter[["COX2"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "COX2")
PMSfilter[["COX1"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "COX1")
PMSfilter[["ATP8"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ATP8")
PMSfilter[["ATP6"]] <- PercentageFeatureSet(PMSfilter, assay = "Spatial", feature = "ATP6")

PMSfilter[["percent.mt"]] <- apply(PMSfilter@meta.data[,5:17], 1, sum)

# p1 <- VlnPlot(PMSfilter, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0.1) + NoLegend()
# plot.name <- paste0("plot/",mysample_name,"_rawcounts_gene_percent.mt-vlnplot.pdf")
# ggsave(filename = plot.name, plot=p1, width = 8, height =6)

#step1: keep protein coding, no chrY genes; use gene name when available, otherwise use gene id
genes.use <- read.table('/Lustre03/data/wangrui/gouyuwei/spatial_transcriptome/protein_gene.txt')
genes.use <- as.vector(genes.use$V1)
PMSfilter_f1<- subset(PMSfilter,features=genes.use)
my_dat[2,1] <- 2      #TODO
my_dat[2,2] <- dim(PMSfilter_f1)[1]    #TODO   
my_dat[2,3] <- dim(PMSfilter_f1)[2]   #TODO


#step 2
min.spots <- 15
num.spots <- rowSums(as.matrix(PMSfilter_f1@assays$Spatial@counts) > 0)
genes.use <- names(num.spots[which(num.spots >= min.spots)])
##mykeepgene <- c(1:nrow(PMSfilter_f1))[rownames(PMSfilter_f1)%in%as.character(genes.use)]
##PMSfilter_f2 <- subset(PMSfilter_f1,features=mykeepgene)
PMSfilter_f2 <- subset(PMSfilter_f1,features=genes.use)
my_dat[3,1] <- 3      #TODO
my_dat[3,2] <- dim(PMSfilter_f2)[1]    #TODO   
my_dat[3,3] <- dim(PMSfilter_f2)[2]   #TODO


#step3
feature_nums <- colSums(as.matrix(PMSfilter_f2@assays$Spatial@counts)>0)
#mycut_feature <- as.numeric(quantile(feature_nums,0.02))
count_nums <- colSums(as.matrix(PMSfilter_f2@assays$Spatial@counts))
#mycut_count<-as.numeric(quantile(count_nums,0.02))
#PMSfilter_f1 <- subset(PMSfilter, subset = nFeature_Spatial>=mycut_feature|nCount_Spatial>=mycut_count& percent.mt < 25)

PMSfilter_f3 <- subset(PMSfilter_f2, subset = nCount_Spatial>=500)
PMSfilter_f3 <- subset(PMSfilter_f3, subset = nFeature_Spatial>=200)
my_dat[4,1] <- 4      #TODO
my_dat[4,2] <- dim(PMSfilter_f3)[1]    #TODO   
my_dat[4,3] <- dim(PMSfilter_f3)[2]   #TODO


#summarize spot number, total gene number, gene number per spot after filter

p1 <- VlnPlot(PMSfilter_f3, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0.1) + NoLegend()
plot.name <- paste(res_path,"p1.pdf",sep = "/")
ggsave(filename = plot.name, plot=p1, width = 8, height =6)

saveRDS(PMSfilter_f3, file = paste(res_path,"data.rds",sep = "/"))     ###TODO

#plot SpatialDimPlot

p2 <- SpatialDimPlot(PMSfilter_f3)
plot.name <- paste(res_path,"p2.pdf",sep = "/")  ####TODO
ggsave(filename = plot.name, plot=p2, width = 8, height =8)


write.table(my_dat,paste(res_path,"count.txt",sep = '/'),quote = F,sep = '\t',row.names = F,col.names = F)

