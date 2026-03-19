
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

	def __init__ (self,path='/Lustre03/data/tangqianzi/Gut_unroll/outputs/'):

		self.path = path
		self.jobfolder = 'runall_09-0stepjobs'
		#self.input = open('/Lustre01/tangqianzi/forLinyu/pwm_files.txt','r')
		self.scriptpath = path

	def generate (self):

                ## alltissues=['BF4','BM4','BM5','TM3','TM4','TF4','TF5']

                os.system('mkdir -p '+self.path+'/'+self.jobfolder+'/')

                ## myoutput=open(self.path+'/'+self.jobfolder+'/runthis_pipeline.sh','w')
	        myoutput1=open(self.path+'/'+self.jobfolder+'/runthis_dump.sh','w')
	        myoutput2=open(self.path+'/'+self.jobfolder+'/alljobs_dump.sh','w')

                mysamples1=['1D','2J_1','2J_2','3I_1','3I_2','4CA_1','4CA_2','4CA_3','5CM','6PC_1','6PC_2','6PC_4','7DC_1']
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
	                       print>>myoutput2,'''export PATH=/Lustre01/tangqianzi/software/anaconda2new/bin/:$PATH
python /Lustre03/data/tangqianzi/Gut_unroll/outputs/09-0.prepare_plot_JEM.py 6 '''+mysample

	        print>>myoutput1,'/usr/bin/perl /Lustre01/tangqianzi/software/scripts/qsub-sgenew.pl --queue tangqianzi --maxjob 200 --lines 2 --jobprefix step9-0 --convert no --resource nodes=1:ppn=1,mem=20g '+self.path+'/'+self.jobfolder+'/alljobs_dump.sh'

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
