
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/

class generate:

	def __init__ (self,path='/Lustre03/data/tangqianzi/scHiC/cell_metabolism_methods/'):

		self.path = path
		self.jobfolder = 'runall_homojobs'
		## self.maxnum = '18'

	def generate (self):

	        #os.system('python '+self.path+'/combine_results.py '+self.mychr+' '+self.path)

	        os.system('mkdir -p '+self.path+'/'+self.jobfolder)
	        #os.system('mkdir -p '+self.outpath)

	        myoutput=open(self.path+'/'+self.jobfolder+'/runthis_pipeline.sh','w')

	        #os.system('python '+self.scriptpath+'/get_loop_union.new.py '+self.mychr)

	        ## samples=['Liver_BF3','Liver_TM3']

	        samples=[]

	        handle=open('/Lustre03/data/tangqianzi/scHiC/cell_metabolism_methods/new_convertTable.txt','r')

	        for line in handle:
	               line=line.rstrip()
	               parts=line.split()
	               #if parts[1]!='X1D':
	               if 1:
	                       samples.append(parts[1])


	        handle.close()


	        ## myconvert={'BF3':'TM3','TM3':'BF3','BF4':'TM4','TM4':'BF4','BM4':'TF4','TF4':'BM4','BM5':'TF5','TF5':'BM5'}

	        print len(samples)

                #======= get simiHICCUPS ===============================

	        myoutput1=open(self.path+'/'+self.jobfolder+'/runthis_test.sh','w')
	        myoutput2=open(self.path+'/'+self.jobfolder+'/alljobs_test.sh','w')

	        for mysample in samples:

                      print>>myoutput2,'''export PATH=/Lustre01/tangqianzi/software/anaconda2new/bin/:$PATH'''
                      print>>myoutput2,'''python '''+self.path+'/get_homoscore.py '+mysample



	        print>>myoutput1,'/usr/bin/perl /Lustre01/tangqianzi/software/scripts/qsub-sgenew.pl --queue tangqianzi --maxjob 100 --lines 6 --jobprefix run_homo --convert no --resource nodes=1:ppn=1,mem=20g '+self.path+'/'+self.jobfolder+'/alljobs_test.sh'

	        myoutput1.close()
	        myoutput2.close()

	        print>>myoutput,'sh '+self.path+'/'+self.jobfolder+'/runthis_test.sh'

		myoutput.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	generate().generate()

	#generate(mychr=pathandfiles[0]).generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
