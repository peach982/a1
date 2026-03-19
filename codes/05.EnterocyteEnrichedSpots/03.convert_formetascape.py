
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
#/Lustre03/data/tangqianzi/Gut_timepoint_bulkRNAseq/ensembl-pig-add-MYH24.gtf

class generate:

	def __init__ (self):

		self.path = '/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/'
		#self.input1 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/limmaresultsfinal/'+time+'_enrichment.'+loci+'.xls','r')
		##self.input1 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/ANOVAlike_analyses/finalresult.genes.xls','r')
		self.input2 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/pig2human.txt','r')
		self.input6 = open('/data/tangyuanling/STsnRNA/forQinglin20251014/gene2types.xls','r')
		self.folder = 'gene_lists_human'

	def generate (self):

	        myloci=['PColon']
	        mytimes=['pre']
	        mydepths=['2']
	        #myloci=['Jejunum','ACaecum','PColon']
	        #mytimes=['pre','0','33','90','180']
	        #mydepths=['2','3','4','5']

	        piggenekeep={}
	        i=0
	        for line in self.input6:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i!=0 and parts[1]=='protein_coding':
	                       piggenekeep[parts[0]]=0

	               i+=1

                i=0
                alldata={}
                for line in self.input2:
                       line=line.rstrip()
                       parts=line.split('\t')
                       if i!=0 and len(parts)>=9:
                            if parts[1]:
                                    pigname=parts[1]
                            else:
                                    pigname=parts[0]

                            if pigname not in piggenekeep:
                                    continue

                            if pigname not in alldata:
                                    alldata[pigname]=[]

                            alldata[pigname].append([int(parts[8]),float(parts[5]),parts[1],parts[3],parts[2]])


                       i+=1

                count=1
                myconvert={}
                for myname in alldata:
                       alldata[myname].sort()
                       alldata[myname].reverse()

                       if not(re.search('ssc-',alldata[myname][0][2])) and not(re.search('MT-',alldata[myname][0][3])) and alldata[myname][0][2]!=alldata[myname][0][3] and alldata[myname][0][2]!='' and alldata[myname][0][3]!='':
                            #print alldata[myname][0][2],alldata[myname][0][3]
                            count+=1


                       myconvert[myname]=alldata[myname][0][4]

                print len(myconvert)

                os.system('mkdir -p '+self.path+'/'+self.folder+'/')

                for mylocus in myloci:

                       for mydepth in mydepths:
                            handleout1=open(self.path+'/'+self.folder+'/'+mylocus+'_depth'+mydepth+'_pseudotime_pos.txt','w')
                            handleout2=open(self.path+'/'+self.folder+'/'+mylocus+'_depth'+mydepth+'_pseudotime_neg.txt','w')
                            handlein=open(self.path+'/'+mylocus+'_depth'+mydepth+'_pseudotime.csv','r')

                            poslist_pig={}
                            neglist_pig={}

                            k=0
                            for line in handlein:
                                line=line.rstrip()
                                parts=line.split(',')
                                if k!=0:
                                        if re.search('pseudotime',line) and float(parts[3])<0.05:
                                                if float(parts[2])>0:
                                                        poslist_pig[parts[0]]=0
                                                else:
                                                        neglist_pig[parts[0]]=0


                                k+=1


                            poslist_human={}
                            neglist_human={}

                            for c in poslist_pig:
                                    if c in myconvert:
                                            poslist_human[myconvert[c]]=0


                            for c in neglist_pig:
                                    if c in myconvert:
                                            neglist_human[myconvert[c]]=0



                            for c in poslist_human:
                                    print>>handleout1,c


                            for c in neglist_human:
                                    print>>handleout2,c


                            ##print poslist_human


                            handlein.close()
                            handleout1.close()
                            handleout2.close()


                       for mytime in mytimes:
                            handleout1=open(self.path+'/'+self.folder+'/'+mylocus+'_time'+mytime+'_pseudotime_pos.txt','w')
                            handleout2=open(self.path+'/'+self.folder+'/'+mylocus+'_time'+mytime+'_pseudotime_neg.txt','w')
                            handlein=open(self.path+'/'+mylocus+'_time'+mytime+'_pseudotime.csv','r')

                            poslist_pig={}
                            neglist_pig={}

                            k=0
                            for line in handlein:
                                line=line.rstrip()
                                parts=line.split(',')
                                if k!=0:
                                        if re.search('pseudotime',line) and float(parts[3])<0.05:
                                                if float(parts[2])>0:
                                                        poslist_pig[parts[0]]=0
                                                else:
                                                        neglist_pig[parts[0]]=0


                                k+=1


                            poslist_human={}
                            neglist_human={}

                            for c in poslist_pig:
                                    if c in myconvert:
                                            poslist_human[myconvert[c]]=0


                            for c in neglist_pig:
                                    if c in myconvert:
                                            neglist_human[myconvert[c]]=0




                            for c in poslist_human:
                                    print>>handleout1,c


                            for c in neglist_human:
                                    print>>handleout2,c


                            handlein.close()
                            handleout1.close()
                            handleout2.close()

		##self.input1.close()
		self.input2.close()
		self.input6.close()

def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	#(options,pathandfiles) = optparser.parse_args()

	generate().generate()

	#generate().generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
