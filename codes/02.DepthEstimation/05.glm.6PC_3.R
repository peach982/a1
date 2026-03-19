setwd('/Lustre03/data/tangqianzi/Gut_unroll/outputs/6PC_3/')

data3<-read.table(file='/Lustre03/data/tangqianzi/Gut_unroll/outputs/pos2spotID.distancefile.xls',sep='\t',header=T)
results<-read.table(file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.lognorm.6.xls',sep='\t',header=T,row.names=1)

#==================================== regression analysis ==========================================
library(splines)
library(qvalue)

y<-matrix(as.vector(results[1,]),4,1)
x<-matrix(as.vector(data3[2:5,1]),4,1)
testdata<-cbind(x,y)
colnames(testdata)<-c("x","y")
testdata<-as.data.frame(testdata)

testdata$x<-as.vector(as.numeric(testdata$x))
testdata$y<-as.vector(as.numeric(testdata$y))

##this is for raw count
p1<-glm(y ~ ns(x, df = 2), data=testdata, family="gaussian")
#p0=glm(y~1,data=testdata, family="gaussian")

#for negative binomial
#anova(p1, test = "Chisq")

#for gaussian
pvalueall<-anova(p1, test = "F")
pvalue<-pvalueall[[6]][2]

get_pvalue<-function(myindex){

y<-matrix(as.vector(results[myindex,]),4,1)
x<-matrix(as.vector(data3[2:5,1]),4,1)
mydata<-cbind(x,y)
colnames(mydata)<-c("x","y")
mydata<-as.data.frame(mydata)

mydata$x<-as.vector(as.numeric(mydata$x))
mydata$y<-as.vector(as.numeric(mydata$y))

p1<-glm(y ~ ns(x, df = 2), data=mydata, family="gaussian")
pvalueall<-anova(p1, test = "F")
pvalue<-pvalueall[[6]][2]
pvalue

}

result_allps<-sapply(1:nrow(results),function(j) get_pvalue(j))
result_allps<-as.matrix(result_allps)
rownames(result_allps)<-rownames(results)
#FDRs<-p.adjust(as.numeric(result_allps[,1]),method ="hommel")
FDRs<-qvalue(as.numeric(result_allps[,1]))
result_allps<-cbind(result_allps,FDRs$qvalue)
colnames(result_allps)<-c("pvalue","qvalue")

finalresult<-cbind(results,result_allps)

write.table(finalresult,file='finalresult_pseudobulk_relativedistance.6.xls',sep='\t',quote=F,row.names=TRUE,col.names=TRUE)
