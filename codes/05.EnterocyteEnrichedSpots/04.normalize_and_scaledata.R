
##source activate /Lustre03/data/ChenZiyu/miniconda3/envs/monocle3
##source activate /Lustre01/tangqianzi/anaconda3forSeurat4/envs/R4.0.2
##export HDF5_USE_FILE_LOCKING='FALSE'

library(Seurat)

setwd('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/')

expression_matrix <- read.table(file='allloci_enterocytes_merged.counts.txt',sep='\t',header=T,row.names=1)
cell_metadata <- read.table(file='allloci_enterocytes_merged.meta.txt',sep='\t',header=T,row.names=1)
#gene_annotation <- read.table(file='/Lustre03/data/tangqianzi/Gut_unroll/outputs/formonocle3noreps/ACaecum/formonocle3_gene_metadata.X.xls',sep='\t',header=T,row.names=1)

my_seurat <- CreateSeuratObject(counts = expression_matrix, meta.data = cell_metadata)
my_seurat <- NormalizeData(object = my_seurat)
my_seurat <- FindVariableFeatures(my_seurat)
my_seurat <- ScaleData(my_seurat,features=rownames(my_seurat))
mymatrix <- as.matrix(my_seurat@assays$RNA[1:dim(my_seurat@assays$RNA)[1],1:dim(my_seurat@assays$RNA)[2]])
colnames(mymatrix) <- colnames(my_seurat@assays$RNA)
rownames(mymatrix) <- rownames(my_seurat@assays$RNA)

write.table(mymatrix,file='ScaledData_allgenes_allenterocytes.txt',sep='\t',quote=F,row.names=T,col.names=T)
