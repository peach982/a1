# -*- coding: utf-8 -*-

import sys
import os
import re
from optparse import OptionParser
import subprocess
import time
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

	def __init__ (self):

	        path='/Lustre03/data/tangqianzi/Gut_unroll/outputs/'

	        #self.input1 = open('/Lustre01/tangqianzi/Gut_project/new/rawcounts/counts_'+sample+'.csv','r')
	        ##self.input1 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/'+sample+'_rawcounts.txt','r')
	        #self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.absolutedistance_finalgroupsforpseudobulk.xls','r')
	        ##self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/spot2group_relative_distance.'+cutoff+'.xls','r')
	        self.input = open(path+'/info.txt','r')
	        self.output1 = open(path+'/ANOVAlike_analyses/00.info.xls','w')
	        self.output2 = open(path+'/ANOVAlike_analyses/00.merged_table.xls','w')
	        self.path = path

	def autoremove_bam (self):

	        allfolders=[]
	        myconvert={}
	        i=0
	        for line in self.input:
	            line=line.rstrip()
	            parts=line.split()
	            if i!=0:
	                   allfolders.append(parts[0])
	                   myconvert[parts[0]]=parts
	            else:
	                   keepnames=parts

	            i+=1

	        allinfo=[]
	        for myfolder in allfolders:
	            for m in range(1,6):
	                allinfo.append(myfolder+'_g'+str(m))


	        print>>self.output1,'rep'+'\t'+'depth'+'\t'+'\t'.join(keepnames)
	        for myrep in allinfo:
	            mysample='_'.join(myrep.split('_')[0:-1])
	            myrestinfo=myconvert[mysample]
	            mydepth=myrep.split('_')[-1].split('g')[1]
	            print>>self.output1,myrep+'\t'+mydepth+'\t'+'\t'.join(myrestinfo)

                allmatrix={}
                for myrep in allinfo:
                    allmatrix[myrep]={}

                i=0
                uniongenes={}
	        for myfolder in allfolders:
	            handle=open(self.path+'/'+myfolder+'/pos2spotID.relativedistance_pseudobulk.PCG.6.mergereps.xls','r')

                    m=0
	            for line in handle:
	                line=line.rstrip()
	                parts=line.split('\t')
	                if m!=0:
	                    uniongenes[parts[0]]=0
	                    for j in range(1,6):
	                        myrepeach=myfolder+'_g'+str(j)
	                        allmatrix[myrepeach][parts[0]]=parts[j]
	                m+=1
	            i+=1


                print>>self.output2,'\t'.join(allinfo)

                for mygene in uniongenes:
                    sig=0
                    for myrep in allinfo:
                        if mygene not in allmatrix[myrep]:
                            sig=1
                    if sig==0:
                        each=[]
                        for myrep in allinfo:
                            each.append(allmatrix[myrep][mygene])
                        print>>self.output2,mygene+'\t'+'\t'.join(each)

	        self.input.close()
	        self.output1.close()
	        self.output2.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	autoremove_bam().autoremove_bam()
	#autoremove_bam().autoremove_bam()


if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
