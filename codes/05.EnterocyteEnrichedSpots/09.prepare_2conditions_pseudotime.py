
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
from scipy.stats import ranksums
import numpy
import math

#export PATH=/Lustre01/tangqianzi/software/anaconda2new/bin/:$PATH
#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
#/Lustre03/data/tangqianzi/Gut_timepoint_bulkRNAseq/ensembl-pig-add-MYH24.gtf

class generate:

	def __init__ (self):


		self.type = 'depth'
		self.info = '5'
		self.segment1 = 'Jejunum'
		self.segment2 = 'ACaecum'
		self.nu = 'water'

		#self.path = path+'/'
		self.input1 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/allloci_enterocytes_pseudotime.txt','r')
		self.input2 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/'+self.nu+'_ME.txt','r')
		self.input3 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/allloci_enterocytes_merged.meta.txt','r')
		self.output = open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/nutrient_MEs/'+self.nu+'_'+self.type+'_'+self.info+'_'+self.segment1+'_'+self.segment2+'.txt','w')
		#self.name = name

		self.path = '/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/nutrient_MEs/'

	def generate (self):

		sample2val={}
		for line in self.input2:
			line=line.rstrip()
			parts=line.split()
			sample2val[parts[0]]=float(parts[1])


		i=0
		allindex={}
		allkeep={}
		sample2info={}
		sample2locus={}
		for line in self.input3:
			line=line.rstrip()
			parts=line.split('\t')
			if i==0:
				for j in range(0,len(parts)):
					allindex[parts[j]]=j+1

			if i!=0:
				myinfo=parts[allindex[self.type]]
				if myinfo==self.info:
					if parts[allindex['locus']]==self.segment1 or parts[allindex['locus']]==self.segment2:
						allkeep[parts[0]]=0

				sample2info[parts[0]]=myinfo
				sample2locus[parts[0]]=parts[allindex['locus']]

			i+=1


		##print allkeep

		sample2pseudotime={}
		for line in self.input1:
			line=line.rstrip()
			parts=line.split()
			sample2pseudotime[parts[0]]=float(parts[1])



		minval=1e6
		maxval=0

		for mysample in allkeep:
			minval=min(minval,sample2pseudotime[mysample])
			maxval=max(maxval,sample2pseudotime[mysample])

		mystep=(maxval-minval)/10.0

		myintervals=[]
		for m in range(0,10):
			start=minval+m*mystep
			end=minval+(m+1)*mystep
			myintervals.append([start,end])


		##print myintervals


		#scipy.stats.ranksums
		samples1=[]
		samples2=[]

		for sample in allkeep:
			if sample2locus[sample]==self.segment1:
				samples1.append(sample)
			if sample2locus[sample]==self.segment2:
				samples2.append(sample)


		print samples1,samples2

		#allresults=[]
		print>>self.output,'group'+'\t'+'shape'+'\t'+'avg_log2FC'+'\t'+'module'
		for m in range(0,10):
			each1=[]
			for mysample in samples1:
				if sample2pseudotime[mysample]>=myintervals[m][0] and sample2pseudotime[mysample]<myintervals[m][1]:
					each1.append(sample2val[mysample])

			each2=[]
			for mysample in samples2:
				if sample2pseudotime[mysample]>=myintervals[m][0] and sample2pseudotime[mysample]<myintervals[m][1]:
					each2.append(sample2val[mysample])


			mypvalue=ranksums(each1,each2)[1]
			if float(mypvalue)<0.05:
				myshape=21
			else:
				myshape=4

			if numpy.mean(each1)!=0 and numpy.mean(each2)!=0:
				mylog2FC=math.log(numpy.mean(each2)/numpy.mean(each1))/math.log(2)
			else:
				mylog2FC=0


			#allresults.append([])
			print>>self.output,str(m+1)+'\t'+str(myshape)+'\t'+str(mylog2FC)+'\t'+self.type+'_'+self.info+'_'+self.segment1+'_'+self.segment2


		self.input1.close()
		self.input2.close()
		self.input3.close()
		self.output.close()



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
