
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
#5% and 0.5

class generate:

	def __init__ (self,sample='1D',cutoff='6'):

		#self.path = path
		self.input = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.relativedistance_pseudobulk.PCG.'+cutoff+'.xls','r')
		self.input2 = open('/Lustre02/data/hic/Gut_project/raw_counts/mito.txt','r')
		self.output = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.relativedistance_pseudobulk.PCG.filtered.'+cutoff+'.xls','w')
		self.output2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.relativedistance_pseudobulk.PCG.filtered.info.'+cutoff+'.xls','w')

	def generate (self):

                mitogenes=[]
	        for line in self.input2:
	               line=line.rstrip()
	               mitogenes.append(line)

	        print>>self.output2,'Cell'+'\t'+'Cell_Cycle'+'\t'+'Locus'+'\t'+'Replicate'

                i=0
                count=0
	        for line in self.input:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i==0:
	                       print>>self.output,line

	                       for j in range(1,len(parts)):
	                               myname=parts[j]
	                               myinfo=parts[j].split('g')[1]
	                               print>>self.output2,myname+'\t'+myname+'\t'+'S'+'\t'+myinfo+'\t'+'1'

	               else:
	                       count_5per=0
	                       allvals=[]
	                       for j in range(1,len(parts)):
	                               if float(parts[j])>0:
	                                       count_5per+=1
	                               allvals.append(float(parts[j]))

	                       mymean=sum(allvals)/float(len(allvals))

	                       if float(count_5per)/float(len(parts)-1)>0.05 and mymean>0.5 and (parts[0] not in mitogenes):
	                               print>>self.output,line
	                               count+=1

	               i+=1

                print count

		self.input.close()
		self.input2.close()
		self.output.close()
		self.output2.close()

def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	#generate().generate()
	generate(cutoff=pathandfiles[0],sample=pathandfiles[1]).generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
