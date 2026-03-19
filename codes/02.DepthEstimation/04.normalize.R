#export R_LIBS=/Lustre01/husilu/R/x86_64-redhat-linux-gnu-library/3.4/:$R_LIBS
#export R_LIBS=/Lustre01/tangqianzi/software/Rlibssilu/:$R_LIBS

library("scater")

setwd('/Lustre03/data/tangqianzi/Gut_unroll/outputs/1D/')
data1<-read.table(file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.6.xls',sep='\t',header=T,row.names=1)
data2<-read.table(file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.info.6.xls',sep='\t',header=T,row.names=1)

example_sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(data1)),
    colData = data2)

example_sce <- normalize(example_sce) # also takes log

#drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) < 3})
#length(drop_genes)-sum(drop_genes==TRUE)
#example_sce <- example_sce[!drop_genes, ]

#example_sce <- runPCA(example_sce)
#example_sce <- runTSNE(example_sce, perplexity=5)

#results<-reducedDim(example_sce, "TSNE")

#write.table(results,file='results_combined_pseudobulk.X.sampe.filtered.xls',sep='\t',quote=F,row.names=T,col.names=F)

results<-as.matrix(example_sce@assays$data[[2]])
colnames(results)<-colnames(example_sce)
rownames(results)<-rownames(example_sce)

write.table(results,file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.lognorm.6.xls',sep='\t',quote=F,row.names=T,col.names=T)

#data3<-read.table(file='pos2spotID.distancefile.xls',sep='\t',header=T)

#==================================== regression analysis ==========================================

library(splines)
y<-matrix(as.vector(results[1,]),10,1)
x<-matrix(as.vector(data3[,1]),10,1)
testdata<-cbind(x,y)
colnames(testdata)<-c("x","y")
testdata<-as.data.frame(testdata)

testdata$x<-as.vector(as.numeric(testdata$x))
testdata$y<-as.vector(as.numeric(testdata$y))

##this is for raw count
p1<-glm(y ~ ns(x, df = 3), data=testdata, family="gaussian")
#p0=glm(y~1,data=testdata, family="gaussian")

#for negative binomial
#anova(p1, test = "Chisq")

#for gaussian
#anova(p1, test = "F")
