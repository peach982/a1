
#source activate /Lustre01/wangrui/anaconda3forSeurat4/envs/R4

args = commandArgs(trailingOnly=TRUE)
samplename<-args[1]

library(STutility)
library(hdf5r)
infoTable<-read.table(paste('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/infoTables/',samplename,'_infoTable.txt',sep=''),sep='\t',header=T)
se <- InputFromTable(infotable = as.data.frame(infoTable),
                      minUMICountsPerGene = 0,
                      minSpotsPerGene = 0,
                      minUMICountsPerSpot = 0,
                      minGenesPerSpot = 0,
                      platform =  "Visium")

se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

celltype<-read.table(paste('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/celltypelabel/',samplename,'_celltype_labels.xls',sep=''),sep='\t',header=T)

se.subset <- SubsetSTData(se, spots = celltype[,1])

orderrownames<-rownames(se.subset@meta.data)

celltype2<-celltype[celltype$cell_id%in%orderrownames,]
rownames(celltype2)<-celltype2$cell_id

celltype_order<-celltype2[orderrownames,2]

se.subset@meta.data$seurat_clusters<-celltype_order

result<-GetSpatNet(se.subset,nNeighbours = 6)
write.table(result[[1]],file=paste('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_GetSpatNet/test_',samplename,'_GetSpatNet.subset.xls',sep=''),sep='\t',quote=F,row.names=T,col.names=T)
write.table(se.subset@meta.data,file=paste('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_metadata/test_',samplename,'_metadata.subset.xls',sep=''),sep='\t',quote=F,row.names=T,col.names=T)

se.subset <- SetIdent(se.subset, value = "seurat_clusters")

#mykeepnames<-c('TA','Paneth','Bcell','Enterocyte','Goblet','EEC','Stem','Myeloid','T_NKcell','Plasma')###ADULT 9
#mykeepnames<-c('Bcell','Enterocyte','Goblet','EEC','Stem','Myeloid','T_NKcell','Plasma')###CHILD 8
#mykeepnames<-c('InnateLymphoid','AdaptiveLymphoid','Enterocyte','Goblet','EEC','Stem','Myeloid','SecretoryProgenitors') ##D0small 8
mykeepnames<-c('InnateLymphoid','AdaptiveLymphoid','Enterocyte','Goblet','EEC','Stem','Myeloid') ##D0big 7

mykeep<-matrix(0,7,7)
rownames(mykeep)<-mykeepnames
colnames(mykeep)<-mykeepnames

for (mid in mykeepnames){

myid<-mid
myid2<-paste('nbs_',myid,sep='')

se2 <- RegionNeighbours(se.subset, id = myid , keep.idents = TRUE, verbose = TRUE)
selfnum<-sum(se2@meta.data[[myid2]][!is.na(se2@meta.data[[myid2]])]==myid)

mykeepnames2<-mykeepnames[ !mykeepnames == myid]

mykeep[myid,myid]<-selfnum
for (c in mykeepnames2){

nonNAnames<-se2@meta.data[[myid2]][!is.na(se2@meta.data[[myid2]])]

nonNAnames2<-c()
for (d in nonNAnames){

nonNAnames2<-c(nonNAnames2,strsplit(d, split = "_")[[1]][1])}

internum<-sum(nonNAnames2==c)
mykeep[myid,c]<-internum

}
}


mykeepall<-vector(mode='list', length=200)
for (i in 1:200){

print (i)
mykeepeach<-matrix(0,7,7)
rownames(mykeepeach)<-mykeepnames
colnames(mykeepeach)<-mykeepnames

mykeepall[[i]]<-mykeepeach

se.subset@meta.data$seurat_clusters<-sample(celltype_order)
se.subset <- SetIdent(se.subset, value = "seurat_clusters")

for (mid in mykeepnames){

myid<-mid
myid2<-paste('nbs_',myid,sep='')

se2 <- RegionNeighbours(se.subset, id = myid , keep.idents = TRUE, verbose = TRUE)
selfnum<-sum(se2@meta.data[[myid2]][!is.na(se2@meta.data[[myid2]])]==myid)

mykeepnames2<-mykeepnames[ !mykeepnames == myid]

mykeepall[[i]][myid,myid]<-selfnum
for (c in mykeepnames2){

nonNAnames<-se2@meta.data[[myid2]][!is.na(se2@meta.data[[myid2]])]

nonNAnames2<-c()
for (d in nonNAnames){

nonNAnames2<-c(nonNAnames2,strsplit(d, split = "_")[[1]][1])}

internum<-sum(nonNAnames2==c)
mykeepall[[i]][myid,c]<-internum

}
}


}


mykeepfinal<-matrix(0,7,7)
rownames(mykeepfinal)<-mykeepnames
colnames(mykeepfinal)<-mykeepnames

for (c in mykeepnames){
for (d in mykeepnames){

formymean<-c()
for (i in 1:200){

formymean<-c(formymean,mykeepall[[i]][c,d])

}

mymean<-mean(formymean)
mySD<-sd(formymean)

mykeepfinal[c,d]<-(mykeep[c,d]-mymean)/mySD

}
}

write.table(mykeepfinal,file=paste('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_heteroscore/test_',samplename,'_heteroscore.subset.xls',sep=''),sep='\t',quote=F,row.names=T,col.names=T)
