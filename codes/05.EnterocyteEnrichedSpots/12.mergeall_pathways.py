
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

		self.output1=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/extract_metascape/merged_pos.xls','w')
		self.output2=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/extract_metascape/merged_neg.xls','w')
		self.output3=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/extract_metascape/merged_both.xls','w')

	def generate (self):

                loci=['Jejunum','ACaecum','PColon']
                mylayers=[]
                mygs=['pos','neg']

                for m in range(2,6):
                        mylayers.append('depth'+str(m))

                mytimes=['pre','0','33','90','180']

                noretain=['PColon_timepre_neg','PColon_timepre_pos','PColon_time0_pos','PColon_time0_neg']

                for mytime in mytimes:
                        mylayers.append('time'+mytime)

                unionterms_pos={}
                unionterms_neg={}
                unionterms_both={}
                alldata_pos={}
                alldata_neg={}
                alldata_both={}

                for mylocus in loci:
                        for mylayer in mylayers:
                                if mylocus+'_'+mylayer+'_pos' not in noretain:
                                        alldata_pos[mylocus+'_'+mylayer+'_pos']={}
                                        handle=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/extract_metascape/'+mylocus+'_'+mylayer+'_pos.txt','r')
                                        m=0
                                        for line in handle:
                                                line=line.rstrip()
                                                parts=line.split('\t')
                                                if m==0:
                                                        m+=1
                                                        continue

                                                if float(parts[-1])<(-1.3):
                                                        unionterms_pos[parts[1]+';'+parts[2]+';'+parts[3]]=0

                                                alldata_pos[mylocus+'_'+mylayer+'_pos'][parts[1]+';'+parts[2]+';'+parts[3]]=0

                                                m+=1

                                        handle.close()


                for mylocus in loci:
                        for mylayer in mylayers:
                                if mylocus+'_'+mylayer+'_neg' not in noretain:
                                        alldata_neg[mylocus+'_'+mylayer+'_neg']={}
                                        handle=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/extract_metascape/'+mylocus+'_'+mylayer+'_neg.txt','r')
                                        m=0
                                        for line in handle:
                                                line=line.rstrip()
                                                parts=line.split('\t')
                                                if m==0:
                                                        m+=1
                                                        continue

                                                if float(parts[-1])<(-1.3):
                                                        unionterms_neg[parts[1]+';'+parts[2]+';'+parts[3]]=0

                                                alldata_neg[mylocus+'_'+mylayer+'_neg'][parts[1]+';'+parts[2]+';'+parts[3]]=0

                                                m+=1

                                        handle.close()


                for mylocus in loci:
                        for mylayer in mylayers:
                                for myg in mygs:
                                        if mylocus+'_'+mylayer+'_'+myg not in noretain:
                                                alldata_both[mylocus+'_'+mylayer+'_'+myg]={}
                                                handle=open('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/extract_metascape/'+mylocus+'_'+mylayer+'_'+myg+'.txt','r')
                                                m=0
                                                for line in handle:
                                                        line=line.rstrip()
                                                        parts=line.split('\t')
                                                        if m==0:
                                                                m+=1
                                                                continue

                                                        if float(parts[-1])<(-1.3):
                                                                unionterms_both[parts[1]+';'+parts[2]+';'+parts[3]]=0

                                                        alldata_both[mylocus+'_'+mylayer+'_'+myg][parts[1]+';'+parts[2]+';'+parts[3]]=0

                                                        m+=1

                                                handle.close()


                printpos=[]
                keeppos=[]
                for mylocus in loci:
                        for mylayer in mylayers:
                                printpos.append(mylocus+'_'+mylayer)
                                keeppos.append(mylocus+'_'+mylayer+'_'+'pos')

                print>>self.output1,'Category'+'\t'+'Term'+'\t'+'Description'+'\t'+'count'+'\t'+'\t'.join(printpos)
                allresults=[]
                for myterm in unionterms_pos:
                        count=0
                        each=[]
                        for mykeep in keeppos:
                                if mykeep not in noretain:
                                      if myterm in alldata_pos[mykeep]:
                                            count+=1
                                            each.append('1')

                                      else:
                                            each.append('0')

                                else:
                                      each.append('na')


                        #print>>self.output1,myterm.split(';')[0]+'\t'+myterm.split(';')[1]+'\t'+myterm.split(';')[2]+'\t'+str(count)+'\t'+'\t'.join(each)
                        myline=myterm.split(';')[0]+'\t'+myterm.split(';')[1]+'\t'+myterm.split(';')[2]+'\t'+str(count)+'\t'+'\t'.join(each)
                        allresults.append([count,myline])


                allresults.sort()
                allresults.reverse()

                for parts in allresults:
                        print>>self.output1,parts[1]


                printneg=[]
                keepneg=[]
                for mylocus in loci:
                        for mylayer in mylayers:
                                printneg.append(mylocus+'_'+mylayer)
                                keepneg.append(mylocus+'_'+mylayer+'_'+'neg')


                print>>self.output2,'Category'+'\t'+'Term'+'\t'+'Description'+'\t'+'count'+'\t'+'\t'.join(printneg)
                allresults=[]
                for myterm in unionterms_neg:
                        count=0
                        each=[]
                        for mykeep in keepneg:
                                if mykeep not in noretain:
                                      if myterm in alldata_neg[mykeep]:
                                            count+=1
                                            each.append('1')

                                      else:
                                            each.append('0')

                                else:
                                      each.append('na')


                        myline=myterm.split(';')[0]+'\t'+myterm.split(';')[1]+'\t'+myterm.split(';')[2]+'\t'+str(count)+'\t'+'\t'.join(each)
                        allresults.append([count,myline])


                allresults.sort()
                allresults.reverse()

                for parts in allresults:
                        print>>self.output2,parts[1]

                printboth=[]
                keepboth=[]
                for myg in mygs:
                        for mylocus in loci:
                                for mylayer in mylayers:
                                        printboth.append(mylocus+'_'+mylayer+'_'+myg)
                                        keepboth.append(mylocus+'_'+mylayer+'_'+myg)


                print>>self.output3,'Category'+'\t'+'Term'+'\t'+'Description'+'\t'+'count'+'\t'+'\t'.join(printboth)
                allresults=[]
                for myterm in unionterms_both:
                        count=0
                        each=[]
                        for mykeep in keepboth:
                                if mykeep not in noretain:
                                      if myterm in alldata_both[mykeep]:
                                            count+=1
                                            each.append('1')

                                      else:
                                            each.append('0')

                                else:
                                      each.append('na')


                        myline=myterm.split(';')[0]+'\t'+myterm.split(';')[1]+'\t'+myterm.split(';')[2]+'\t'+str(count)+'\t'+'\t'.join(each)
                        allresults.append([count,myline])



                allresults.sort()
                allresults.reverse()

                for parts in allresults:
                        print>>self.output3,parts[1]

		self.output1.close()
		self.output2.close()
		self.output3.close()


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
