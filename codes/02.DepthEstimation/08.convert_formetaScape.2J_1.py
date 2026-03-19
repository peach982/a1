
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
#/Lustre03/data/tangqianzi/Gut_timepoint_bulkRNAseq/ensembl-pig-add-MYH24.gtf

class generate:

	def __init__ (self,path='/Lustre03/data/tangqianzi/Gut_myproject_new/limmaresultsfinal/',name='2J_1'):

		#self.path = path+'/forDAVID/'
		#self.input1 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/limmaresultsfinal/'+time+'_enrichment.'+loci+'.xls','r')
		self.input1 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+name+'/J_one_step_countresults.cluster.cut0.5.txt','r')
		self.input2 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/limmaresultsfinal/pig2human.txt','r')
		self.input6 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/limmaresultsfinal/gene2types.xls','r')
		self.output1 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/formetaScape/'+name+'/increase_formetaScape.txt','w')
		self.output2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/formetaScape/'+name+'/decrease_formetaScape.txt','w')
		self.output3 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/formetaScape/'+name+'/midhigh_formetaScape.txt','w')
		self.output4 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/formetaScape/'+name+'/increasedrop_formetaScape.txt','w')

	def generate (self):

	        piggenekeep={}
	        i=0
	        for line in self.input6:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i!=0 and parts[1]=='protein_coding':
	                       piggenekeep[parts[0]]=0

	               i+=1


                i=0
                increase=[]
                decrease=[]
                midhigh=[]
                midlow=[]
                #middle={}
                for line in self.input1:
                       line=line.rstrip()
                       parts=line.split('\t')
                       if i!=0:
                            if parts[1]=='1' or parts[1]=='4':
                                    increase.append(parts[0])
                            if parts[1]=='7' or parts[1]=='8':
                                    decrease.append(parts[0])
                            if parts[1]=='3' or parts[1]=='5':
                                    midhigh.append(parts[0])
                            if parts[1]=='2':
                                    midlow.append(parts[0])


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
                            print alldata[myname][0][2],alldata[myname][0][3]
                            count+=1


                       myconvert[myname]=alldata[myname][0][4]

                print count


                keepnames1={}
                keepnames2={}
                keepnames3={}
                keepnames4={}
                for myname in increase:
                       if myname in myconvert:
                            keepnames1[myconvert[myname]]=0
                       else:
                            #print myname
                            #print "wrong"
                            pass

                for myname in decrease:
                       if myname in myconvert:
                            keepnames2[myconvert[myname]]=0
                       else:
                            #print myname
                            #print "wrong"
                            pass

                for myname in midhigh:
                       if myname in myconvert:
                            keepnames3[myconvert[myname]]=0
                       else:
                            #print myname
                            #print "wrong"
                            pass


                for myname in midlow:
                       if myname in myconvert:
                            keepnames4[myconvert[myname]]=0
                       else:
                            #print myname
                            #print "wrong"
                            pass

                for c in keepnames1:
                       print>>self.output1,c

                for c in keepnames2:
                       print>>self.output2,c

                for c in keepnames3:
                       print>>self.output3,c

                for c in keepnames4:
                       print>>self.output4,c

                print len(keepnames1)

		self.input1.close()
		self.input2.close()
		self.input6.close()
		self.output1.close()
		self.output2.close()
		self.output3.close()
		self.output4.close()


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
