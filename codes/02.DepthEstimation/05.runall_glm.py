
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import gzip

#cd runall_fimorestjobs
#nohup sh runthis_dump.sh &

class generate:

	def __init__ (self,path='/Lustre03/data/tangqianzi/Gut_unroll/outputs/',cutoff='6'):

		self.path = path
		self.jobfolder = 'runall_step05jobs'
		#self.input = open('/Lustre01/tangqianzi/forLinyu/pwm_files.txt','r')
		#self.scriptpath = path
		self.cutoff = cutoff

	def generate (self):

                ## alltissues=['BF4','BM4','BM5','TM3','TM4','TF4','TF5']

                os.system('mkdir -p '+self.path+'/'+self.jobfolder+'/')

                ## myoutput=open(self.path+'/'+self.jobfolder+'/runthis_pipeline.sh','w')
	        myoutput1=open(self.path+'/'+self.jobfolder+'/runthis_dump.sh','w')
	        myoutput2=open(self.path+'/'+self.jobfolder+'/alljobs_dump.sh','w')

                mysamples1=['1D','2J_1','2J_2','3I_1','3I_2','4CA_1','4CA_2','4CA_3','5CM','6PC_1','6PC_2','6PC_3','6PC_4','7DC_1']
                mysamples2=['P0_2_M_1','P0_3_K_3','P0_4_J_2','P33_4_J_1','P33_5_K_1','P33_5_M_1','P90_2_J_2','P90_2_K_3','P90_2_M_1','P90_4_M_1']

                mysamples=[]
                for c in mysamples1:
                        mysamples.append(c)



                ##mysamples=[]
                for c in mysamples2:
                        mysamples.append(c)

                print len(mysamples)

                ## for mytissue in alltissues:
                if 1:
                      ##os.system('mkdir -p '+self.path+'/'+mytissue+'/')
	              #for mymotif in allmotifs:
	              for mysample in mysamples:
	                       #os.system('mkdir -p /Lustre01/wangrui/HiC_GWAS/noreplace/'+mymotif+'/')
	                       #os.system('mkdir -p /Lustre01/tangqianzi/simulations_NingNC/chicken_genome/'+mychr+'/')

	                       handle=open(self.path+'/'+self.jobfolder+'/'+mysample+'.R','w')

	                       print>>handle,'''setwd("/Lustre03/data/tangqianzi/Gut_unroll/outputs/'''+mysample+'''/")
data3<-read.table(file='/Lustre03/data/tangqianzi/Gut_unroll/outputs/pos2spotID.distancefile.xls',sep='\t',header=T)
results<-read.table(file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.lognorm.'''+self.cutoff+'''.xls',sep='\t',header=T,row.names=1)

#==================================== regression analysis ==========================================
library(splines)
library(qvalue)

y<-matrix(as.vector(results[1,]),5,1)
x<-matrix(as.vector(data3[,1]),5,1)
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
pvalueall<-anova(p1, test = "F")
pvalue<-pvalueall[[6]][2]

get_pvalue<-function(myindex){

y<-matrix(as.vector(results[myindex,]),5,1)
x<-matrix(as.vector(data3[,1]),5,1)
mydata<-cbind(x,y)
colnames(mydata)<-c("x","y")
mydata<-as.data.frame(mydata)

mydata$x<-as.vector(as.numeric(mydata$x))
mydata$y<-as.vector(as.numeric(mydata$y))

p1<-glm(y ~ ns(x, df = 3), data=mydata, family="gaussian")
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

write.table(finalresult,file='finalresult_pseudobulk_relativedistance.'''+self.cutoff+'''.xls',sep='\t',quote=F,row.names=TRUE,col.names=TRUE)
'''


	                       handle.close()

	                       print>>myoutput2,'''export R_LIBS=/Lustre01/husilu/R/x86_64-redhat-linux-gnu-library/3.4/:$R_LIBS
export R_LIBS=/Lustre01/tangqianzi/software/Rlibssilu/:$R_LIBS
/usr/bin/Rscript '''+self.path+'/'+self.jobfolder+'/'+mysample+'.R'

	        print>>myoutput1,'/usr/bin/perl /Lustre01/tangqianzi/software/scripts/qsub-sgenew.pl --queue tangqianzi --maxjob 200 --lines 3 --jobprefix step5 --convert no --resource nodes=1:ppn=1,mem=30g '+self.path+'/'+self.jobfolder+'/alljobs_dump.sh'

	        myoutput1.close()
	        myoutput2.close()

	        ## print>>myoutput,'sh '+self.path+'/'+self.jobfolder+'/runthis_dump.sh'

	        ## myoutput.close()
	        #self.input.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	generate().generate()



if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
