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
#randomly sample 20 spots for each distance
#dupplicated ids allowed here!!!

class autoremove_bam:

	def __init__ (self,sample='1D',cutoff='6'):

	        #self.input1 = open('/Lustre01/tangqianzi/Gut_project/new/rawcounts/counts_'+sample+'.csv','r')
	        self.input1 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/'+sample+'_rawcounts.txt','r')
	        #self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.absolutedistance_finalgroupsforpseudobulk.xls','r')
	        self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/spot2group_relative_distance.'+cutoff+'.xls','r')
	        self.input3 = open('/Lustre02/data/hic/Gut_project/raw_counts/pseudobulk.X/gene2types.xls','r')
	        self.output = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.relativedistance_pseudobulk.PCG.'+cutoff+'.xls','w')
                self.samplesize = '100'

	def autoremove_bam (self):

	        keepgenes={}
	        for line in self.input3:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if parts[1]=='protein_coding':
	                       keepgenes[parts[0]]=1


                i=0
                group2ids={}
                allids={}
	        for line in self.input2:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i!=0:
	                       if parts[0] not in group2ids:
	                               group2ids[parts[0]]=[]

	                       group2ids[parts[0]].append(parts[1])
	                       allids[parts[0]]=0

	               i+=1


                allids_list=[]
                for c in allids:
                        allids_list.append(int(c))

                allids_list.sort()

                i=0
                for line in self.input1:
                        line=line.rstrip()
                        parts=line.split('\t')
                        if i==0:
                                name2index={}
                                for j in range(0,len(parts)):
                                        myname=parts[j]
                                        name2index[myname]=j+1

                                group2index={}
                                printnames=[]
                                for c in allids_list:
                                        eachindex=[]
                                        for myid in group2ids[str(c)]:
                                                eachindex.append(name2index[myid])

                                        group2index[str(c)]=eachindex

                                        printnames.append('g'+str(c))

                                print>>self.output,'geneid'+'\t'+'\t'.join(printnames)

                                group2index2={}
                                for c in allids_list:
                                        each=[]
                                        for k in range(0,int(self.samplesize)):
                                                each.append(random.sample(group2index[str(c)],1)[0])

                                        group2index2[str(c)]=each


                        else:

                                printvals=[]
                                for c in allids_list:
                                        eachsum=0
                                        for j in group2index2[str(c)]:
                                                eachsum+=int(parts[j])

                                        printvals.append(str(eachsum))

                                newgenename=parts[0]

                                if newgenename in keepgenes:
                                        print>>self.output,newgenename+'\t'+'\t'.join(printvals)

                        i+=1

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

	autoremove_bam(cutoff=pathandfiles[0],sample=pathandfiles[1]).autoremove_bam()


if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
