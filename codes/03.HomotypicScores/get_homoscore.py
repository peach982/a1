
import sys
import os
import re
import numpy as np
from optparse import OptionParser
import subprocess
import time
import random

#min 2 circles
#max 4 circles: chr 14
#shuffled = random.sample(a_list, len(a_list))

class autoremove_bam:

	def __init__ (self,name='X2J_1'):

		self.input1 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_GetSpatNet/test_'+name+'_GetSpatNet.subset.xls','r')
		self.input2 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_metadata/test_'+name+'_metadata.subset.xls','r')
		self.output = open('/data/tangyuanling/STsnRNA/forQinglin20251014/01.homoscoresTQZ/outputs_homoscore/test_'+name+'_homoscore.xls','w')

		self.pertimes = '200'

	def autoremove_bam (self):

	        celltypes=['TA','Paneth','Bcell','Enterocyte','Goblet','EEC','Stem','Myeloid','T_NKcell','Plasma']
#mykeepnames<-c('TA','Paneth','Bcell','Enterocyte','Goblet','EEC','Stem','Myeloid','T_NKcell','Plasma')###ADULT 9
#mykeepnames<-c('Bcell','Enterocyte','Goblet','EEC','Stem','Myeloid','T_NKcell','Plasma')###CHILD 8
#mykeepnames<-c('InnateLymphoid','AdaptiveLymphoid','Enterocyte','Goblet','EEC','Stem','Myeloid','SecretoryProgenitors') ##D0small 8
#mykeepnames<-c('InnateLymphoid','AdaptiveLymphoid','Enterocyte','Goblet','EEC','Stem','Myeloid') ##D0big 7

	        cellid2type={}
	        keepids=[]
	        keeptypes=[]
	        i=0
	        for line in self.input2:
	                line=line.rstrip()
	                parts=line.split('\t')
	                if i!=0:
	                       cellid2type[parts[0]]=parts[6]
	                       keepids.append(parts[0])
	                       keeptypes.append(parts[6])

	                i+=1


                allkeeptypes={}
                for i in range(0,int(self.pertimes)):
                        allkeeptypes[str(i)]=[]

	        for m in range(0,int(self.pertimes)):
	                keeptypes_each=random.sample(keeptypes, len(keeptypes))
	                allkeeptypes[str(m)]=keeptypes_each

                alldata=[]
                allcellid2type={}
                allresults={}
                for i in range(0,int(self.pertimes)):
                        allcellid2type[str(i)]={}
                        for m in range(0,len(keepids)):
                            myid=keepids[m]
                            mytype=allkeeptypes[str(i)][m]

                            allcellid2type[str(i)][myid]=mytype








                i=0
                for line in self.input1:
                        line=line.rstrip()
                        parts=line.split('\t')
                        alldata.append(parts)
                        if i!=0:
                                #myfromcell=cellid2type[parts[1]]
                                #mytocell=cellid2type[parts[2]]
                                #if parts[1] not in allresults:
                                allresults[parts[1]]=[]
                                #allresults[parts[1]].append(parts[2])
                                #if myfromcell==mytocell:
                                        #allresults[myfromcell][parts[1]]+=1


                        i+=1


                i=0
                for parts in alldata:
                        if i!=0:
                                allresults[parts[1]].append(parts[2])
                                #if myfromcell==mytocell:
                                        #allresults[myfromcell][parts[1]]+=1


                        i+=1



                obscore={}
                expscore={}

                for celltype in celltypes:
                    obscore[str(celltype)]=[]

                for m in range(0,int(self.pertimes)):
                    expscore[str(m)]={}
                    for celltype in celltypes:
                        expscore[str(m)][celltype]=[]

                for myfromid in allresults:
                        myfromcell=cellid2type[myfromid]
                        countsame=0
                        for c in allresults[myfromid]:
                            mytocell=cellid2type[c]
                            if mytocell==myfromcell:
                                countsame+=1


                        if myfromcell in celltypes:
                                obscore[str(myfromcell)].append(countsame)

                        for m in range(0,int(self.pertimes)):
                            myfromcell=allcellid2type[str(m)][myfromid]
                            countsame=0
                            for c in allresults[myfromid]:
                                mytocell=allcellid2type[str(m)][c]
                                if mytocell==myfromcell:
                                    countsame+=1

                            if myfromcell in celltypes:
                                    expscore[str(m)][myfromcell].append(countsame)




                for celltype in celltypes:

                    if len(obscore[str(celltype)])!=0:
                            myob=sum(obscore[str(celltype)])/float(len(obscore[str(celltype)]))
                            myexpall=[]
                            for m in range(0,int(self.pertimes)):
                                print float(len(expscore[str(m)][celltype]))
                                print celltype
                                myexpeach=sum(expscore[str(m)][celltype])/float(len(expscore[str(m)][celltype]))
                                myexpall.append(myexpeach)

                            myexp=sum(myexpall)/float(len(myexpall))
                            print>>self.output,str(celltype)+'\t'+str(myob)+'\t'+str(myexp)+'\t'+str(myob-myexp)

                    else:
                            print>>self.output,str(celltype)+'\t'+'NA'+'\t'+'NA'+'\t'+'NA'



		self.input1.close()
		self.input2.close()
		self.output.close()

def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	autoremove_bam(name=pathandfiles[0]).autoremove_bam()


if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
