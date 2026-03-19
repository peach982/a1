##source activate /Lustre03/data/ChenZiyu/miniconda3/envs/monocle3

library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)
library(data.table)

setwd('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/')

expression_matrix <- read.table(file='allloci_enterocytes_merged.counts.txt',sep='\t',header=T,row.names=1)
cell_metadata <- read.table(file='allloci_enterocytes_merged.meta.txt',sep='\t',header=T,row.names=1)
gene_annotation <- read.table(file='/data/tangyuanling/STsnRNA/forQinglin20251014/formonocle3noreps/Jejunum/formonocle3_gene_metadata.X.xls',sep='\t',header=T,row.names=1)

# Make the CDS object
cds <- new_cell_data_set(as.matrix(expression_matrix),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "batch")

#cds <- align_cds(cds, residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

cds <- reduce_dimension(cds)
##plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "depth")
cds <- cluster_cells(cds)
cds <- learn_graph(cds) ######,use_partition=FALSE)

get_earliest_principal_node <- function(cds, time_bin="1"){
  cell_ids <- which(colData(cds)[, "depth"] == time_bin & colData(cds)[, "time"] == -45)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

traj.coord<- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
peudoresult<-as.matrix(traj.coord)

pdf("depth.enterocytes.100.pdf")
plot_cells(cds,
           color_cells_by = "depth",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, graph_label_size=4, cell_size = 1, trajectory_graph_color = "darkred", trajectory_graph_segment_size = 0.5)
dev.off()


pdf("time.enterocytes.100.pdf")
plot_cells(cds,
           color_cells_by = "time",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, graph_label_size=4, cell_size = 1, trajectory_graph_color = "darkred", trajectory_graph_segment_size = 0.5)
dev.off()


pdf("locus.enterocytes.100.pdf")
plot_cells(cds,
           color_cells_by = "locus",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, graph_label_size=4, cell_size = 1, trajectory_graph_color = "darkred", trajectory_graph_segment_size = 0.5)
dev.off()

pdf("pseudotime.enterocytes.100.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=4, cell_size = 1, trajectory_graph_color = "black", trajectory_graph_segment_size = 0.5)
dev.off()

write.table(peudoresult,file='allloci_enterocytes_pseudotime.txt',sep='\t',quote=F,row.names=T,col.names=F)

#============ !!! left for analyses after regression !!! ===========================
##compare pre70 and the rest of time points: from layer 1 to layer 5
##vioplot: 5 colors, darkred, darkblue, darkgreen, orange, purple

mylist<-c(-45,0,33,90,180)
mylist2<-c('pre','0','33','90','180')
names(mylist2)<-mylist
mylocus<-'ACaecum'

for (k in 1:5){
mytime<-mylist[k]
mytime2<-mylist2[k]

cell_ids <- which(colData(cds)[, "locus"] == mylocus & colData(cds)[, "time"] == mytime)
cds_subset <- cds[,cell_ids]
gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime", expression_family="negbinomial")
fit_coefs <- coefficient_table(gene_fits)
test_result<-fit_coefs %>% select(gene_short_name, term, estimate, q_value)
fwrite(test_result,file=paste(mylocus,'_time',mytime2,'_pseudotime.csv',sep=''),quote=F,row.names=F,col.names=T)

}


for (k in 2:5){
mydepth<-k

cell_ids <- which(colData(cds)[, "locus"] == mylocus & colData(cds)[, "depth"] == mydepth)
cds_subset <- cds[,cell_ids]
gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime", expression_family="negbinomial")
fit_coefs <- coefficient_table(gene_fits)
test_result<-fit_coefs %>% select(gene_short_name, term, estimate, q_value)
fwrite(test_result,file=paste(mylocus,'_depth',mydepth,'_pseudotime.csv',sep=''),quote=F,row.names=F,col.names=T)

}


#cell_ids <- which(colData(cds)[, "locus"] == 'Jejunum' & colData(cds)[, "time"] == -45)
#cds_subset <- cds[,cell_ids]
#gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime", expression_family="negbinomial")
#fit_coefs <- coefficient_table(gene_fits)
#test_result<-fit_coefs %>% select(gene_short_name, term, estimate, q_value)
#fwrite(test_result,file='/Lustre03/data/tangqianzi/test.csv',quote=F,row.names=F,col.names=T)


##gene_fits1 <- fit_models(cds_subset, model_formula_str = "~pseudotime+batch", expression_family="negbinomial")
##fit_coefs1 <- coefficient_table(gene_fits1)
##mytable1<-as.matrix(fit_coefs1)
##write.table(mytable1,file='Jejunum_pre_pseudotime.time_batch.txt',sep='\t',quote=F,row.names=F,col.names=T)

##checkres<-compare_models(gene_fits1, gene_fits) %>% select(gene_short_name, q_value)
