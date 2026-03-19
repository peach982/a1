
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

	def __init__ (self):

		#self.path = path+'/'
		self.input2 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/JEM_transport.new.txt','r')
		self.input1 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/ScaledData_allgenes_allenterocytes.humansymbols.txt','r')
		#self.name = name
		self.path = '/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/'

	def generate (self):

	        allinfo={}

	        for line in self.input2:
	               line=line.rstrip()
	               parts=line.split('\t')
	               mynutrient='_'.join(parts[1].split(' '))
	               mygene=parts[0]

	               if mynutrient not in allinfo:
	                       allinfo[mynutrient]=[]

	               allinfo[mynutrient].append(mygene)



	        alldata={}
	        k=0
	        for line in self.input1:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if k==0:
	                       headline='\t'.join(parts[2:-1])
	               else:

	                       alldata[parts[0]]='\t'.join(parts[2:-1])
	               k+=1

	        for mynutrient in allinfo:
	               handle=open(self.path+'/enterocyte_matrix_'+mynutrient+'.txt','w')

	               print>>handle,headline

	               for mygene in allinfo[mynutrient]:
	                       if mygene in alldata:
	                               print>>handle,mygene+'\t'+alldata[mygene]

	               handle.close()
	               ##handle


		self.input1.close()
		self.input2.close()


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
