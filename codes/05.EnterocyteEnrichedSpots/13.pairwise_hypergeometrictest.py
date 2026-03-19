
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
		self.output=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/pairwise_hypergeometric.xls','w')

	def generate (self):

                loci=['Jejunum','ACaecum','PColon']
                mylayers=[]
                mygs=['pos','neg']

                for m in range(2,6):
                        mylayers.append('depth'+str(m))

                mytimes=['pre','0','33','90','180']

                noretain=['PColon_timepre_neg','PColon_timepre_pos','PColon_time0_pos','PColon_time0_neg','PColon_time33_neg']

                for mytime in mytimes:
                        mylayers.append('time'+mytime)


                allgenesets={}
                allbackgenes={}
                for mylocus in loci:
                        for mylayer in mylayers:
                                allgenesets[mylocus+'_'+mylayer+'_pos']={}
                                allgenesets[mylocus+'_'+mylayer+'_neg']={}
                                handle=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/'+mylocus+'_'+mylayer+'_pseudotime.csv','r')

                                m=0
                                for line in handle:
                                        line=line.rstrip()
                                        parts=line.split(',')
                                        if m!=0:
                                                if re.search('pseudotime',line) and float(parts[-1])<0.05:
                                                        if float(parts[-2])<0:
                                                                allgenesets[mylocus+'_'+mylayer+'_neg'][parts[0]]=0

                                                        else:
                                                                allgenesets[mylocus+'_'+mylayer+'_pos'][parts[0]]=0
                                                allbackgenes[parts[0]]=0

                                        m+=1

                                handle.close()

                allranks=[]
                for myg in mygs:
                        for mylocus in loci:
                                for mylayer in mylayers:
                                        allranks.append(mylocus+'_'+mylayer+'_'+myg)


                print>>self.output,'info'+'\t'+'\t'.join(allranks)
                for m in range(0,len(allranks)):
                        each=[]
                        for n in range(0,len(allranks)):
                                if m==n:
                                        pvalue='0'
                                else:
                                        geneset1=allgenesets[allranks[m]]
                                        geneset2=allgenesets[allranks[n]]

                                        a1=0
                                        for gene in geneset1:
                                                if gene in geneset2:
                                                        a1+=1

                                        a1=a1-1

                                        a2=len(geneset1)
                                        a3=0

                                        for gene in allbackgenes:
                                                if gene not in geneset1:
                                                        a3+=1

                                        a4=len(geneset2)

                                        pvalue=str(1-hypergeom.cdf(a1, a2+a3, a4, a2, loc=0))

                                each.append(pvalue)

                        print>>self.output,allranks[m]+'\t'+'\t'.join(each)

#phyper(10-1, 16, 13, 11, lower.tail=F)
#1-hypergeom.cdf(10-1, 16+13, 11, 16, loc=0)

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
