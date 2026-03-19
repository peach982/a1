
argv<-commandArgs(trailingOnly=T)
input_nu<-argv[1]

setwd('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/')
folder_path <- paste('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/nutrient_MEs/pdfs/',input_nu,'/',sep='')
dir.create(folder_path)

library(ggplot2)

data1<-read.table(file='/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/allloci_enterocytes_merged.meta.txt',sep='\t',header=T,row.names=1)
data2<-read.table(file=paste(input_nu,'_ME','.txt',sep=''),sep='\t',header=F,row.names=1)
data3<-read.table(file='allloci_enterocytes_pseudotime.txt',sep='\t',header=F,row.names=1)

mydata<-cbind(data2[,1],data3[,1])
colnames(mydata)<-c('ME','pseudotime')

rownames(mydata)<-rownames(data2)

#Jejunm depth 2~5
#Jejunm time -45,0,33,90,180

mylist<-c(-45,0,33,90,180)
mylist2<-c('pre','0','33','90','180')
myloci<-c('Jejunum','ACaecum','PColon')

for (index1 in 1:3){

mylocus<-myloci[index1]

for (index2 in 2:5){
mydepth<-index2
keepcells<-rownames(data1[data1$locus==mylocus&data1$depth==mydepth,])
mydata2<-mydata[rownames(mydata)%in%keepcells,]
mydata2<-as.data.frame(mydata2)
print(dim(mydata2))

p <- ggplot(mydata2, aes(pseudotime, ME)) +
  geom_point()

# loess method: local regression fitting
p <- p + geom_smooth(method = "loess")

ggsave(paste("/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/nutrient_MEs/pdfs/",input_nu,"/",input_nu,"_",mylocus,"_depth",mydepth,"_MEandpseudotime.pdf",sep=""), plot = p)

}

for (index3 in 1:5){
mytime<-mylist[index3]
mytimename<-mylist2[index3]

keepcells<-rownames(data1[data1$locus==mylocus&data1$time==mytime,])
mydata2<-mydata[rownames(mydata)%in%keepcells,]
mydata2<-as.data.frame(mydata2)
print(dim(mydata2))

p <- ggplot(mydata2, aes(pseudotime, ME)) +
  geom_point()

# loess method: local regression fitting
p <- p + geom_smooth(method = "loess")

ggsave(paste("/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/nutrient_MEs/pdfs/",input_nu,"/",input_nu,"_",mylocus,"_time",mytimename,"_MEandpseudotime.pdf",sep=""), plot = p)


}
}
