
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
		self.jobfolder = 'runall_step04jobs'
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

	                       print>>handle,'''library("scater")

setwd("/Lustre03/data/tangqianzi/Gut_unroll/outputs/'''+mysample+'''/")
data1<-read.table(file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.'''+self.cutoff+'''.xls',sep='\t',header=T,row.names=1)
data2<-read.table(file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.info.'''+self.cutoff+'''.xls',sep='\t',header=T,row.names=1)

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

write.table(results,file='pos2spotID.relativedistance_pseudobulk.PCG.filtered.lognorm.'''+self.cutoff+'''.xls',sep='\t',quote=F,row.names=T,col.names=T)
'''


	                       handle.close()

	                       print>>myoutput2,'''export R_LIBS=/Lustre01/husilu/R/x86_64-redhat-linux-gnu-library/3.4/:$R_LIBS
export R_LIBS=/Lustre01/tangqianzi/software/Rlibssilu/:$R_LIBS
/usr/bin/Rscript '''+self.path+'/'+self.jobfolder+'/'+mysample+'.R'

	        print>>myoutput1,'/usr/bin/perl /Lustre01/tangqianzi/software/scripts/qsub-sgenew.pl --queue tangqianzi --maxjob 200 --lines 3 --jobprefix step4 --convert no --resource nodes=1:ppn=1,mem=30g '+self.path+'/'+self.jobfolder+'/alljobs_dump.sh'

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
