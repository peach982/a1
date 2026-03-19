
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import gzip

#cd runall_fimorestjobs
#nohup sh runthis_dump.sh &

#2,5: lower
#3,4,6: higher
#1: mid higher
#7,9: mid lower

class generate:

	def __init__ (self,path='/Lustre03/data/tangqianzi/Gut_unroll/outputs/',cutoff='6'):

		self.path = path
		self.jobfolder = 'runall_step07_6PC_3jobs'
		#self.input = open('/Lustre01/tangqianzi/forLinyu/pwm_files.txt','r')
		self.scriptpath = path
		self.cutoff = cutoff

	def generate (self):

                ## alltissues=['BF4','BM4','BM5','TM3','TM4','TF4','TF5']

                os.system('mkdir -p '+self.path+'/'+self.jobfolder+'/')


                ## myoutput=open(self.path+'/'+self.jobfolder+'/runthis_pipeline.sh','w')
	        myoutput1=open(self.path+'/'+self.jobfolder+'/runthis_dump.sh','w')
	        myoutput2=open(self.path+'/'+self.jobfolder+'/alljobs_dump.sh','w')

                #mysamples=['1D','2J_1','2J_2','3I_1','3I_2','4CA_1','4CA_2','4CA_3','5CM','6PC_1','6PC_2','6PC_3','6PC_4','7DC_1']

                mysamples=['6PC_3']

                ## for mytissue in alltissues:
                if 1:
                      ##os.system('mkdir -p '+self.path+'/'+mytissue+'/')
	              #for mymotif in allmotifs:
	              for mysample in mysamples:

	                       os.system('mkdir -p '+'/Lustre03/data/tangqianzi/Gut_unroll/outputs/formetaScape/'+mysample+'/')
	                       #os.system('mkdir -p /Lustre01/wangrui/HiC_GWAS/noreplace/'+mymotif+'/')
	                       #os.system('mkdir -p /Lustre01/tangqianzi/simulations_NingNC/chicken_genome/'+mychr+'/')

	                       handle=open(self.path+'/'+self.jobfolder+'/'+mysample+'.R','w')

	                       print>>handle,'''library(edgeR)
library(DESeq2)
library(maSigPro)
library(MASS)

setwd("/Lustre03/data/tangqianzi/Gut_unroll/outputs/'''+mysample+'''/")

data_group <- read.table(file="pos2spotID.distancefile.masigpro.'''+self.cutoff+'''.xls",sep="\t",header=T,row.names=1)
data <- read.table(file="finalresult_pseudobulk_relativedistance.masigpro.'''+self.cutoff+'''.xls",sep="\t",header=T,row.names=1)

#dds <- DESeqDataSetFromMatrix(countData = data,
#                             colData = data_group,
#                             design= ~ Time)

#dds <- DESeq(dds)
#dds <- estimateSizeFactors(dds)
#norm.data<-counts(dds,normalized=TRUE)

norm.data<-data

norm.data.new<-norm.data[ rowSums(norm.data)!=0, ]

#data_design<-read.table(file="design_file.txt",sep="\t",header=T,row.names=1)

#data_design <- data_group
#data_design <- as.matrix(data_design)
#rownames(data_design) <- rownames(data_group)
#colnames(data_design) <- c("Time","Replicate")
data_design<-as.data.frame(data_group)

##data_design$Replicate<-c(rep(c(1,2,3,4),5))

##data_design$Tissue <- c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
d <- make.design.matrix(data_design)

NBp <- p.vector(norm.data, d, counts=TRUE, min.obs=4)
NBt <- T.fit(NBp)

#get<-get.siggenes(NBt, vars="groups", rsq = 0.5)
get2<-get.siggenes(NBt, vars="all", rsq = 0.5)

pdf("/Lustre03/data/tangqianzi/Gut_unroll/outputs/formetaScape/'''+mysample+'''/pseudobulkPG_Result.0.5.count.pdf")
testresult2<-see.genes(get2$sig.genes, show.fit = T, dis = d$dis, newX11=FALSE)
dev.off()

newcut<-as.matrix(testresult2$cut)

write.table(get2$sig.genes$sig.pvalues,file="J_one_step_countresults.cut0.5.txt",sep="\t",quote=F,row.names=T,col.names=T)
write.table(newcut,file="J_one_step_countresults.cluster.cut0.5.txt",sep="\t",quote=F,row.names=T,col.names=T)
'''


	                       handle.close()

	                       print>>myoutput2,'''export R_LIBS=/Lustre01/husilu/R/x86_64-redhat-linux-gnu-library/3.4/:$R_LIBS
export R_LIBS=/Lustre01/tangqianzi/software/Rlibssilu/:$R_LIBS
Rscript '''+self.path+'/'+self.jobfolder+'/'+mysample+'.R'

	        print>>myoutput1,'/usr/bin/perl /Lustre01/tangqianzi/software/scripts/qsub-sgenew.pl --queue tangqianzi --maxjob 200 --lines 3 --jobprefix step7 --convert no --resource nodes=1:ppn=1,mem=30g '+self.path+'/'+self.jobfolder+'/alljobs_dump.sh'

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
