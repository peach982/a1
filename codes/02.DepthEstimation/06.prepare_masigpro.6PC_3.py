
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/

##P0,P180,P32,P90,Pe70
##J,K,M

class generate:

	def __init__ (self,sample='6PC_3',cutoff='6'):

	        self.path = '/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/'

		self.input1 = open(self.path+'/finalresult_pseudobulk_relativedistance.'+cutoff+'.xls','r')
		self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/pos2spotID.distancefile.xls','r')

		self.output1 = open(self.path+'/finalresult_pseudobulk_relativedistance.masigpro.'+cutoff+'.xls','w')
		self.output2 = open(self.path+'/pos2spotID.distancefile.masigpro.'+cutoff+'.xls','w')

	def generate (self):

                i=0
                for line in self.input1:
                        line=line.rstrip()
                        parts=line.split('\t')
                        if i==0:
                               print>>self.output1,'\t'.join(parts[:4])

                        else:
                               #if float(parts[-2])<0.05:
                               if float(parts[-2])<0.15:
                                    count0=0
                                    for c in parts[1:5]:
                                            if float(c)==0:
                                                    count0+=1

                                    if count0<2:
                                            print>>self.output1,'\t'.join(parts[:5])


                        i+=1

                i=0
                for line in self.input2:
                        line=line.rstrip()
                        parts=line.split('\t')
                        if i==0:
                             print>>self.output2,'Time'+'\t'+'Replicate'+'\t'+'Group'

                        elif i>=2:

                             print>>self.output2,'g'+str(i)+'\t'+str(parts[0])+'\t'+'1'+'\t'+'1'

                        i+=1


		self.input1.close()
		self.input2.close()
		self.output1.close()
		self.output2.close()



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
