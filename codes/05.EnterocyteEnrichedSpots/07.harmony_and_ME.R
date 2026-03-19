
library(harmony)
library(tibble)

setwd('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/')
term<-c('sugar','Nucleotide','water','vitamin_absorption','metal_ion','lipid','inorganic_solutes','bile_salt','organic_solutes','amino_acid')
for (i in term){
expression_matrix <- read.table(file=paste0('enterocyte_matrix_',i,'.txt'),sep='\t',header=T,row.names=1)
cell_metadata <- read.table(file='allloci_enterocytes_merged.meta.txt',sep='\t',header=T,row.names=1)

cell_metadata$cell_id<-rownames(cell_metadata)
newtib <- cell_metadata %>%
    as_tibble(rownames = "cell_id")

harmony_embeddings <- HarmonyMatrix(
    t(expression_matrix), newtib, vars_use='batch', do_pca = TRUE, verbose=TRUE, npcs = 3, max_iter = 50, return_object=TRUE
)

mykeep<-matrix(t(harmony_embeddings$Z_corr[1,]),ncol=1)
rownames(mykeep)<-rownames(cell_metadata)

write.table(mykeep,file=paste0(i,'_ME.txt'),sep='\t',quote=F,row.names=T,col.names=F)
}