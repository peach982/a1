
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import numpy
import math
from scipy.stats import hypergeom

#export PATH=/Lustre01/tangqianzi/software/anaconda2new/bin/:$PATH
#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
#/Lustre03/data/tangqianzi/Gut_timepoint_bulkRNAseq/ensembl-pig-add-MYH24.gtf

class generate:

	def __init__ (self):

                self.path='/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/'

		self.input1=open(self.path+'/allloci_enterocytes_merged.meta.txt','r')
		self.input2=open(self.path+'/allloci_enterocytes_pseudotime.txt','r')



	def generate (self):


                output=open(self.path+'/plot_pseudotime_vioplot.R','w')

                #first by loci, then by time, then by depth
                myloci=['Jejunum','ACaecum','PColon']
                mytimes=['-45','0','33','90','180']
                mydepths=['1','2','3','4','5']

                #different scale 5 depth
                #5 times: darkred, darkgreen, darkblue, purple, orange
                m=0
                allsamples={}
                for line in self.input1:
                        line=line.rstrip()
                        parts=line.split('\t')
                        if m!=0:
                                myinfo=parts[-2]+'_'+parts[-1]+'_'+parts[1]
                                if myinfo not in allsamples:
                                        allsamples[myinfo]=[]

                                allsamples[myinfo].append(parts[0])


                        m+=1


                sample2time={}
                for line in self.input2:
                        line=line.rstrip()
                        parts=line.split('\t')
                        sample2time[parts[0]]=float(parts[1])


                labels=[]
                ##colors=[]
                names=[]
                alldata={}
                count=1

                allcolors=['darkred','darkgreen','darkblue','purple','orange']

                print>>output,'allCol<-c()'

                for k in range(0,3):
                        for m in range(0,5):
                                print>>output,'Col<-colorRampPalette(c("white","'+allcolors[m]+'"))(6)'
                                print>>output,'allCol<-c(allCol,Col[2:6])'


                for mylocus in myloci:
                        ##countsig=0
                        for mytime in mytimes:
                                ##myoricolor=allcolors[countsig]
                                for mydepth in mydepths:
                                        alldata[mylocus+'_'+mytime+'_'+mydepth]=[]
                                        if mylocus+'_'+mytime+'_'+mydepth in allsamples:
                                                for c in allsamples[mylocus+'_'+mytime+'_'+mydepth]:
                                                      alldata[mylocus+'_'+mytime+'_'+mydepth].append(str(sample2time[c]))


                                        labels.append('x'+str(count))
                                        print>>output,'x'+str(count)+'<-c('+','.join(alldata[mylocus+'_'+mytime+'_'+mydepth])+')'

                                        names.append('"'+mylocus+'_'+mytime+'_'+mydepth+'"')

                                        count+=1



                                ##countsig+=1

                print>>output,'''pdf("'''+self.path+'''/plot_pseudotime_vioplot.pdf",height=5,width=15)
library(vioplot)
par(mar = c(15,5,3,3))
mylist<-c('''+','.join(labels)+''')
mymatrix<-list('''+','.join(labels)+''')
mylabels<-c('''+','.join(names)+''')
mypos<-c(1:(3*5*5))
mypos2<-c(1:(3*5*5))
mycolors<-allCol
myborders<-allCol
mycolMed<-c(rep("gray",(3*5*5)))

plot(c(),c(),xlim=c(0,(3*5*5)+1),ylim=range(mylist),type="n",ylab="Pseudotime",xlab="",axes=FALSE)

#axis(side=1,at=mypos,labels=mylabels)

axis(side=1,at=c(1:(3*5*5)), labels=F, tick=F)

text(c(1:(3*5*5)), -2.3, srt = 45, adj = 1, labels = mylabels, xpd = TRUE, cex=0.7)

axis(side=2)

for (i in 1:(3*5*5)){

vioplot(mymatrix[[i]],at=mypos2[i],col=mycolors[i],add=TRUE,border=myborders[i],colMed=mycolMed[i])


}

##legend("bottomright",c("Low-altitude","High-altitude"), fill=c("white","gray"), bg="white", text.col = "black", cex=1, bty="n")
dev.off()
'''


		output.close()
		self.input1.close()
		self.input2.close()

		os.system('Rscript '+self.path+'/plot_pseudotime_vioplot.R')



def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	#generate().generate()
	generate().generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
