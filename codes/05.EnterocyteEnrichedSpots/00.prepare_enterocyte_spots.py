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

#random.sample(a,3) default: without replacement/repeated samples

class autoremove_bam:

	def __init__ (self,sample='Pre70_Jejunum',cutoff='6'):

	        path='/data/tangyuanling/STsnRNA/forQinglin20251014/'

	        self.input3 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/gene2types.xls','r')
                self.input1 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/convert_names.new.txt','r')

	        ##self.output = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.relativedistance_pseudobulk.PCG.'+cutoff+'.mergereps.xls','w')

	        self.output1 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/'+sample+'/enterocyte_spots.expression.sampled.xls','w')
	        self.output2 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/'+sample+'/enterocyte_spots.meta.sampled.xls','w')

                self.samplesize = '100'
                self.path = path
                self.sample = sample
                #self.cutratio = '0.5'

	def autoremove_bam (self):

	        gene2type={}
	        for line in self.input3:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if parts[1]=='protein_coding':
	                       gene2type[parts[0]]=1

                keepnames=[]
	        for line in self.input1:
	               line=line.rstrip()
	               parts=line.split(' ')
	               if parts[2]==self.sample:
	                       keepnames.append(parts[0])

                sample2groupall={}
                #sample2batch={}
                count=0
                allrankedsamplenames=[]
                keepgenes=[]
                alldata={}
	        for myname in keepnames:
                        print(myname)
                        handle1=open('/data/tangyuanling/STsnRNA/forQinglin20251014/rawcounts/counts_'+myname+'.csv','r')
                        handle2=open(self.path+'/'+myname+'/spot2group_relative_distance.6.xls','r')
                        handle3=open('/data/tangyuanling/STsnRNA/forQinglin20251014/deconvolve/'+myname+'.labeltransfer3.xls','r')

                        keep_onecelltype={}
                        k=0
                        for line in handle3:
                                line=line.rstrip()
                                parts=line.split('\t')
                                if k==0:
                                    for j in range(1,len(parts)):
                                        if parts[j]=='predicted.id':
                                                id_index=j+1

                                        #if parts[j]=='prediction.score.max':
                                                #score_index=j+1


                                else:
                                        #myscore=float(parts[score_index])
                                        myid=parts[id_index]
                                        if myid=='Enterocyte':
                                        #if myscore>=float(self.cutratio) and myid=='Enterocyte':
                                                ##print "hi"
                                                myname2=myname+'-'+parts[0]
                                                myname2='X'+re.sub('-','_',myname2)
                                                keep_onecelltype[myname2]=0
                                                ##print myname2,"1"

                                k+=1

                        k=0
                        for line in handle2:
                               line=line.rstrip()
                               parts=line.split('\t')
                               if k!=0:
                                    myname2=myname+'-'+parts[1]
                                    myname2='X'+re.sub('-','_',myname2)
                                    sample2groupall[myname2]=parts[0]

                                    #print myname2

                                    #sample2batch[myname2]=str(count+1)

                               k+=1

                        i=0
                        for line in handle1:
                                line=line.rstrip()
                                parts=line.split(',')
                                if i==0:
                                        #headnames=[]
                                        keepindex=[]
                                        for j in range(1,len(parts)):
                                                newname=parts[j].split('"')[1]
                                                newname=re.sub('\.','-',newname)
                                                ##newname=myname+'-'+parts[j].split('"')[1]
                                                newname=myname+'-'+newname
                                                newname='X'+re.sub('-','_',newname)

                                                ##print newname,"2"

                                                if (newname in sample2groupall) and (newname in keep_onecelltype):

                                                        ##print "hi"

                                                        #headnames.append(newname)
                                                        allrankedsamplenames.append(newname)
                                                        keepindex.append(j)

                                        #if count==0:
                                                #headline='\t'.join(headnames)
                                        #else:
                                                #headline=headline+'\t'+'\t'.join(headnames)

                                else:
                                        genename=parts[0].split('"')[1]

                                        neweach=[]

                                        for j in keepindex:
                                                neweach.append(parts[j])

                                        if genename in gene2type:

                                                if gene2type[genename]==1:
                                                        if count==0:
                                                            keepgenes.append(genename)

                                                if gene2type[genename]==1:
                                                        if count==0:
                                                                alldata[genename]=neweach
                                                        else:
                                                                finaleach=[]
                                                                for c in alldata[genename]:
                                                                        finaleach.append(c)

                                                                for c  in neweach:
                                                                        finaleach.append(c)

                                                                ##alldata[genename]=alldata[genename]+'\t'+'\t'.join(neweach)
                                                                alldata[genename]=finaleach


                                i+=1


                        handle1.close()
                        handle2.close()

                        count+=1


                allgroups2index={}
                selectedgroups2index={}
                for m in range(1,6):
                        allgroups2index[str(m)]=[]


                print len(allrankedsamplenames),"test"

                for j in range(0,len(allrankedsamplenames)):
                        mysample=allrankedsamplenames[j]
                        myg=sample2groupall[mysample]
                        allgroups2index[str(myg)].append(j)


                ##print allgroups2index['1']

                printnames_depth=[]
                printnames_samples=[]
                printnames_batches=[]
                for m in range(1,6):
                        selectedgroups2index[str(m)]=[]
                        #eachindex=[]
                        #for k in range(0,int(self.samplesize)):
                        #        each.append(random.sample(allgroups2index[str(m)],1)[0])

                        if int(self.samplesize)<=len(allgroups2index[str(m)]):
                                 eachindex=random.sample(allgroups2index[str(m)],int(self.samplesize))
                        else:
                                 eachindex=allgroups2index[str(m)]


                        print len(eachindex)

                        for k in range(0,len(eachindex)):
                                 printnames_depth.append(str(m))

                        for k in eachindex:
                                 printnames_samples.append(allrankedsamplenames[k])
                                 eachname='_'.join(allrankedsamplenames[k].split('_')[:-2])
                                 printnames_batches.append(eachname)

                        selectedgroups2index[str(m)]=eachindex


                ##print>>self.output,'geneid'+'\t'+'\t'.join(printnames)

                print>>self.output1,'\t'.join(printnames_samples)
                print>>self.output2,'depth'+'\t'+'batch'

                for gene in keepgenes:
                        myvals=alldata[gene]

                        printvals=[]
                        for m in range(1,6):
                                eachsum=[]
                                for j in selectedgroups2index[str(m)]:
                                        eachsum.append(str(myvals[j]))

                                printvals.extend(eachsum)

                        print>>self.output1,gene+'\t'+'\t'.join(printvals)


                for k in range(0,len(printnames_samples)):
                        print>>self.output2,printnames_samples[k]+'\t'+printnames_depth[k]+'\t'+printnames_batches[k]

	        self.input1.close()
	        #self.input2.close()
	        self.input3.close()
	        self.output1.close()
	        self.output2.close()

def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	autoremove_bam(sample=pathandfiles[0]).autoremove_bam()



if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
