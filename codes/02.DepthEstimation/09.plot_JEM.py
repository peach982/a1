
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#export PATH=/Lustre01/tangqianzi/software/anaconda2new/bin/:$PATH
#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
#/Lustre03/data/tangqianzi/Gut_timepoint_bulkRNAseq/ensembl-pig-add-MYH24.gtf

class generate:

	def __init__ (self,path='/Lustre03/data/tangqianzi/Gut_unroll/outputs/',name='1D',cutoff='6'):

		self.path = path+'/'
		#self.input1 = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/pipeline20220812/'+name+'/finalresult_pseudobulk_relativedistance.6.humanids.xls','r')
		self.input1 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+name+'/finalresult_pseudobulk_relativedistance.'+cutoff+'.humanids.xls','r')
		self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/JEM_transport.new.txt','r')
		#self.output1 = open(path+'/forvioplot.'+name+'.JEM.txt','w')
		#self.output2 = open(path+'/forvioplot.'+name+'.JEM.R','w')
		self.name = name

	def generate (self):

                types2genes={}
	        for line in self.input2:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if parts[1] not in types2genes:
	                       types2genes[parts[1]]=[]

	               types2genes[parts[1]].append(parts[0])


	        myranks=['g1','g2','g3','g4','g5']

	        mydatacounts={'g1':1,'g2':1,'g3':1,'g4':1,'g5':1}


	        i=0

	        alldata={}
	        for line in self.input1:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i==0:
	                       ##mynames=parts[:-2]

                               keepnames=[]
	                       for myrank in myranks:
	                               keepnames.append(str(myrank))

	                       keepindex=[]
	                       for myname in keepnames:
	                               for j in range(0,len(parts)):
	                                       if myname==parts[j]:
	                                               keepindex.append(j)



	                       headline='\t'.join(keepnames)

                               #print keepindex
                               #print mydatacounts

	               else:

                               each=[]
                               each.append(parts[0])
                               for j in keepindex:
                                    each.append(str(parts[j]))

                               #contain genename in the each line
                               alldata[parts[0]]='\t'.join(each)

	               i+=1


	        os.system('mkdir -p '+'/Lustre03/data/tangqianzi/Gut_unroll/outputs/JEM/'+self.name+'/')

                allcount=0
	        for mytype in types2genes:
	               finalname='_'.join(mytype.split(' '))
	               handle1=open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/JEM/'+self.name+'/'+finalname+'.txt','w')

	               print>>handle1,headline

	               for mygene in types2genes[mytype]:
	                       if mygene in alldata:
	                               print>>handle1,alldata[mygene]

	               handle1.close()


	               handle2=open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/JEM/'+self.name+'/'+finalname+'.R','w')

	               print>>handle2,'''library(vioplot)
library(scales)

##args = commandArgs(trailingOnly=TRUE)

##myname<-args[1]

#myname<-"amino_acid"

setwd("/Lustre03/data/tangqianzi/Gut_unroll/outputs/JEM/'''+self.name+'''/")

data<-read.table(file=paste("'''+finalname+'''",".txt",sep=""),sep="\t",header=T,row.names=1)

datanew<-c()

for (i in 1:nrow(data)){

myvec<-(as.numeric(data[i,])-mean(as.numeric(data[i,])))/sd(as.numeric(data[i,]))
datanew<-rbind(datanew,myvec)

}

rownames(datanew)<-rownames(data)
colnames(datanew)<-colnames(data)

allcolors1<-c(rep("#CCFF99",'''+str(mydatacounts['g1'])+'''),rep("#99CCCC",'''+str(mydatacounts['g2'])+'''),rep("#66CCFF",'''+str(mydatacounts['g3'])+'''),rep("#990066",'''+str(mydatacounts['g4'])+'''),rep("pink",'''+str(mydatacounts['g5'])+'''))

setwd("/Lustre03/data/tangqianzi/Gut_unroll/outputs/JEM/'''+self.name+'''/")
pdf(paste("'''+finalname+'''",".vioplot.pdf",sep=""),height=5,width=5)
par(mar=c(8,6,5,2))
mylist<-c()
for (i in 1:ncol(datanew)){
mylist<-c(mylist,as.vector(as.numeric(datanew[,i])))
}

plot(c(),c(),xlim=c(0,length(colnames(datanew))+1),ylim=range(mylist),type="n",ylab="Normalized expression levels",xlab="",axes=FALSE,las=2)

newnames<-c()
for (i in 1:length(colnames(datanew))){

#newname<-paste(strsplit(colnames(datanew)[i],"_")[[1]][1],strsplit(colnames(datanew)[i],"_")[[1]][2],sep='_')
newname<-colnames(datanew)[i]
print (newname)
newnames<-c(newnames,newname)

}

axis(side=1,at=c(1:length(colnames(datanew))), labels=newnames, tick=F, las=2)

axis(side=2,las=2)

for (i in 1:length(colnames(datanew))){

vioplot(datanew[,i],at=i,col=alpha(allcolors1[i],0.7),add=TRUE,border=allcolors1[i],colMed="black",lwd=2)

}

medianall<-c()
for (i in 1:ncol(datanew)){
medianall<-c(medianall,round(median(datanew[,i]),1))

}

mtext(c("median=",medianall),side=3,line=1,col="black",at=c(-0.5,1:length(colnames(datanew))),cex=0.8)
dev.off()
'''

                       handle2.close()

                       os.system('export R_LIBS=/Lustre01/husilu/R/x86_64-redhat-linux-gnu-library/3.4/:$R_LIBS')
                       os.system('export R_LIBS=/Lustre01/tangqianzi/software/Rlibssilu/:$R_LIBS')
                       os.system('/usr/bin/Rscript '+'/Lustre03/data/tangqianzi/Gut_unroll/outputs/JEM/'+self.name+'/'+finalname+'.R')

                       allcount+=1

                       #if allcount>=1:
                                #break


		self.input1.close()
		self.input2.close()
		#self.output1.close()
		#self.output2.close()



def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	#generate().generate()

	generate(cutoff=pathandfiles[0],name=pathandfiles[1]).generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
