# -*- coding: utf-8 -*-

import sys
import os
import re
import numpy as np
from optparse import OptionParser
import subprocess
import time
from scipy import stats
import numpy as np
import math
import random

#XX, X_Y
#y – y1={ (y2 – y1)/ (x2 – x1)}(x – x1)
#export PATH=/Lustre01/tangqianzi/software/anaconda2new/bin/:$PATH
#Muscle_TM4.5kb.raw.PEI.xls.keep_PEIs.pick.FDR6_Dis25k.cut_0.2OE.keepbyFrequency
#just test SV and enhancers
#just test SV and promoters
#stats.ttest_ind()[1]
#startid: boundary spots
#endid: in_spots

class autoremove_bam:

	def __init__ (self,pathin='/Lustre03/data/tangqianzi/Gut_unroll/outputs/',cutoff='6',sample='P90_2_K_3'):

	        self.input = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.test.xls','r')
	        #self.input2 = open('/Users/hsl/Desktop/Gut_unroll/unroll2021208/1D_UMAP.new.xls','r')
	        self.input2 = open('/Lustre03/data/wangrui/gouyuwei/spatial_transcriptome/remove_csv/'+sample+'_XYZ_remove.csv','r')
	        self.output1 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.relative_distance.summary.'+cutoff+'.xls','w')
	        self.output2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/spot2group_relative_distance.'+cutoff+'.xls','w')

	        self.cutoff = cutoff

	def autoremove_bam (self):

	        i=0
	        keepids={}
	        for line in self.input2:
	               line=line.rstrip()
	               parts=line.split(',')
	               if i!=0:
	                       if parts[1]=='X':
	                               keepids[parts[0]]=0

	               i+=1


                i=0
                #alldata=[]
                id2dist={}
                bound2spot={}
                alldists0={}
	        for line in self.input:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i!=0 and (parts[1] in keepids):
	               #if i!=0:
	                       #alldata.append(parts)
                               id2dist[parts[1]]=parts[2]
                               #alldists[parts[2]]=0
                               if parts[0] not in bound2spot:
                                    bound2spot[parts[0]]=[]

                               bound2spot[parts[0]].append(parts[1])
                               alldists0[parts[2]]=0

	               i+=1


                allresults={}
                for c in alldists0:
                        if float(c)>=int(self.cutoff):
                                allresults[c]={}
                                for m in range(1,6):
                                        allresults[c][str(m)]=0

                alldists_list=[]
                for c in alldists0:
                        alldists_list.append(float(c))

                alldists_list.sort()

                print>>self.output2,'group_id'+'\t'+'inner_spotid'+'\t'+'base_spotid'
	        for mybound in bound2spot:
	               alldists=[]
	               for c in bound2spot[mybound]:
	                       mydist=float(id2dist[c])
	                       alldists.append(mydist)

	               maxval=max(alldists)
	               #minval=min(alldists)
	               minval=0

	               mystep=(maxval-minval)/5.0

	               allcuts_start=[]
	               allcuts_end=[]

	               for m in range(0,5):
	                       allcuts_start.append(m*mystep+minval)
	                       allcuts_end.append((m+1)*mystep+minval)



	               #if maxval>=5:
	               if maxval>=int(self.cutoff):
	                       for c in bound2spot[mybound]:
	                               mydist=float(id2dist[c])
	                               for m in range(0,4):
	                                       if mydist>=allcuts_start[m] and mydist<allcuts_end[m]:
	                                               mygroup=m+1

	                               if mydist>=allcuts_start[4] and mydist<=allcuts_end[4]:
	                                       mygroup=4+1

	                               allresults[str(maxval)][str(mygroup)]+=1

	                               print>>self.output2,str(mygroup)+'\t'+c+'\t'+mybound


                print>>self.output1,'max_distance'+'\t'+'g1'+'\t'+'g2'+'\t'+'g3'+'\t'+'g4'+'\t'+'g5'

                allsum={}

                for m in range(1,6):
                        allsum[str(m)]=0

                for c in alldists_list:
                        if c>=int(self.cutoff):
                                each=[]
                                for m in range(1,6):
                                        each.append(str(allresults[str(c)][str(m)]))
                                        allsum[str(m)]+=allresults[str(c)][str(m)]

                                print>>self.output1,str(c)+'\t'+'\t'.join(each)


                each=[]
                for m in range(1,6):
                     each.append(str(allsum[str(m)]))

                print>>self.output1,'summation'+'\t'+'\t'.join(each)

                #self.output2.close()
	        self.input.close()
	        self.input2.close()
	        self.output1.close()
	        self.output2.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	autoremove_bam(cutoff=pathandfiles[0],sample=pathandfiles[1]).autoremove_bam()


if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
