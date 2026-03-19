library(Seurat)
setwd("/data/tangyuanling/STsnRNA/forQinglin20251014")
sample<-c("1D","2J_1","2J_2","3I_1","3I_2","4CA_1","4CA_2","4CA_3","5CM","6PC_1","6PC_2","6PC_3","7DC_1","P0_2_M_1", 'P0_3_K_3',"P0_4_J_2",'P33_4_J_1','P33_5_K_1','P33_5_M_1', 'P90_2_J_2','P90_2_K_3','P90_2_M_1','P90_4_M_1',"Pe_1_M_1", "Pe_2_J_1","Pe_1_K_1")
#############predicted.id
for (m in sample) {

  file_path <- paste0('/data/liuqinglin/10.7/spatial/data_spatial/rds/', m, "-filter+SCT+PCA.rds")

  ST <- readRDS(file_path)
  ST <- subset(ST, new == 'X')
df<-read.table(paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,"_PANDA_prop.csv"),sep=',',header=T)
n_col <- ncol(df)
n_row <- nrow(df)
rownames(df)<-df[,1]
df=df[,-1]  
df[, 'predicted.id']<-0  

for (i in 1:nrow(df)) {
  max_value <- max(df[i,-ncol(df)])
  row_data <- df[i, -ncol(df), drop = FALSE]
  
  # 找出与max列值相等（在容差内）的列名
  equal_cols <- colnames(row_data)[row_data[1,]==max_value]
  df[i, 'predicted.id']<-equal_cols
  # 替换最后一列的值
}

df<-as.data.frame(df)
meta<-ST@meta.data
meta$X<-rownames(meta)
head(meta)
head(df)
df$X<-rownames(meta)
meta <- merge(meta, df[, c("X", "predicted.id")], by = "X", all.x = TRUE)
head(meta)
rownames(meta)<-meta[,'x']
write.table(meta,paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,".labeltransfer.xls"),sep='\t',quote=F,row.names=T)
}
#############predicted.id+Enterocyte比例
sample<-c("1D","2J_1","2J_2","3I_1","3I_2","4CA_1","4CA_2","4CA_3","5CM","6PC_1","6PC_2","6PC_3","7DC_1","P0_2_M_1", 'P0_3_K_3',"P0_4_J_2",'P33_4_J_1','P33_5_K_1','P33_5_M_1', 'P90_2_J_2','P90_2_K_3','P90_2_M_1','P90_4_M_1',"Pe_1_M_1", "Pe_2_J_1","Pe_1_K_1")
for (m in sample) {
  file_path <- paste0('/data/liuqinglin/10.7/spatial/data_spatial/rds/', m, "-filter+SCT+PCA.rds")

  ST <- readRDS(file_path)
  ST <- subset(ST, new == 'X')
df<-read.table(paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,"_PANDA_prop.csv"),sep=',',header=T)
n_col <- ncol(df)
n_row <- nrow(df)
rownames(df)<-df[,1]
df=df[,-1]  
df[, 'predicted.id']<-0  

for (i in 1:nrow(df)) {
  max_value <- max(df[i,-ncol(df)])
  row_data <- df[i, -ncol(df), drop = FALSE]
  
  # 找出与max列值相等（在容差内）的列名
  equal_cols <- colnames(row_data)[row_data[1,]==max_value]
  df[i, 'predicted.id']<-equal_cols
  # 替换最后一列的值
}

df<-as.data.frame(df)
meta<-ST@meta.data
meta$X<-rownames(meta)
head(meta)
head(df)
df$X<-rownames(meta)
meta <- merge(meta, df[, c('X',"predicted.id",'Enterocyte')], by = "X", all.x = TRUE)
head(meta)
rownames(meta)<-meta[,'X']
write.table(meta,paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,".labeltransfer2.xls"),sep='\t',quote=F,row.names=T)
}
####查看哪些样本层数差别大
sample<-c("1D","2J_1","2J_2","3I_1","3I_2","4CA_1","4CA_2","4CA_3","5CM","6PC_1","6PC_2","6PC_3","7DC_1","P0_2_M_1", 'P0_3_K_3',"P0_4_J_2",'P33_4_J_1','P33_5_K_1','P33_5_M_1', 'P90_2_J_2','P90_2_K_3','P90_2_M_1','P90_4_M_1',"Pe_1_M_1", "Pe_2_J_1","Pe_1_K_1")
for (m in sample) {
data<-read.table(paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,".labeltransfer2.xls"),sep='\t')
print(m)
print(dim(data))
a<-read.table(paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/",m,"/spot2group_relative_distance.6.xls"),sep='\t',header=T)
print(dim(a))
}
########
#[1] "3I_1"
#[1] 1077   27
#[1] 200   3
#[1] "P0_4_J_2"
#[1] 903  26
#[1] 111   3

#######正常再提取前20%的enterocyte   ##P33_4_J_1是本身Enterocyte就少
sample<-c("1D","2J_1","2J_2","3I_2","4CA_1","4CA_2","4CA_3","5CM","6PC_1","6PC_2","6PC_3","7DC_1","P0_2_M_1", 'P0_3_K_3','P33_4_J_1','P33_5_K_1','P33_5_M_1', 'P90_2_J_2','P90_2_K_3','P90_2_M_1','P90_4_M_1',"Pe_1_M_1", "Pe_2_J_1","Pe_1_K_1")
for (m in sample) {
data<-read.table(paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,".labeltransfer2.xls"),sep='\t')
print(m)
threshold <- quantile(data$Enterocyte, 0.8)  # 获取80%分位数
data2<- data[data$Enterocyte >= threshold, ]
#print(table(data2$predicted.id))
write.table(data2,paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,".labeltransfer3.xls"),sep='\t',quote=F,row.names=T)
b<-read.table(paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/",m,"/spot2group_relative_distance.6.xls"),sep='\t',header=T)
print(length(intersect(data2$X,b$inner_spotid)))
print(dim(data2))
}
########反常样取全部
sample<-c("3I_1","P0_4_J_2")
for (m in sample) {
data<-read.table(paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,".labeltransfer2.xls"),sep='\t')
write.table(data,paste0("/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/",m,".labeltransfer3.xls"),sep='\t',quote=F,row.names=T)
}
