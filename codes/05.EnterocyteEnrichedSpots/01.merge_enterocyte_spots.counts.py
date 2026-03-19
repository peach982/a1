
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import gzip

#cd runall_fimorestjobs
#nohup sh runthis_dump.sh &

class generate:

	def __init__ (self,path='/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/'):

		self.path = path
		#self.jobfolder = 'runall_prepare_enterocyte_spotsjobs'
		self.input = open('/data/tangyuanling/STsnRNA/forQinglin20251014/convert_names.new.2.txt','r')
		self.output1 = open(path+'/allloci_enterocytes_merged.counts.txt','w')
		self.output2 = open(path+'/allloci_enterocytes_merged.meta.txt','w')
		self.oripath = '/data/tangyuanling/STsnRNA/forQinglin20251014/'

	def generate (self):


	        print>>self.output2,'depth'+'\t'+'batch'+'\t'+'locus'+'\t'+'time'

	        myconvert_locus={}
	        myconvert_time={}
	        for line in self.input:
	               line=line.rstrip()
	               parts=line.split(' ')
	               myconvert_locus[parts[2]]=parts[4]
	               myconvert_time[parts[2]]=parts[3]

                mysamples=['Pre70_Jejunum','Pre70_Caecum','Pre70_Colon','P0_Caecum','P0_Jejunum','P0_Colon','P33_Caecum','P33_Jejunum','P33_Colon','P90_Caecum','P90_Jejunum','P90_Colon','ACaecum','Jejunum','PColon']

                alldata={}
                allspots=[]
                allgenes={}
                allgenes_order=[]
                count=0
                for mysample in mysamples:
                        alldata[mysample]={}
                        handle1=open(self.oripath+'/'+mysample+'/enterocyte_spots.expression.sampled.xls','r')
                        handle2=open(self.oripath+'/'+mysample+'/enterocyte_spots.meta.sampled.xls','r')

                        i=0
                        for line in handle1:
                                line=line.rstrip()
                                parts=line.split('\t')
                                if i==0:
                                        for c in parts:
                                                allspots.append(c)


                                else:
                                        alldata[mysample][parts[0]]=parts[1:]

                                        #if count==0:
                                        allgenes[parts[0]]=0

                                        if count==0:
                                                allgenes_order.append(parts[0])


                                i+=1



                        i=0
                        for line in handle2:
                                line=line.rstrip()
                                parts=line.split('\t')
                                if i!=0:
                                        print>>self.output2,line+'\t'+myconvert_locus[mysample]+'\t'+myconvert_time[mysample]

                                i+=1



                        handle1.close()
                        handle2.close()

                        count+=1

	        #myoutput1.close()
	        #myoutput2.close()

	        ## print>>myoutput,'sh '+self.path+'/'+self.jobfolder+'/runthis_dump.sh'

	        print>>self.output1,'\t'.join(allspots)

	        for mygene in allgenes_order:

                       sig=1
	               for mysample in mysamples:
                            if mygene not in alldata[mysample]:
                                sig=0

                       if sig==1:

                           each=[]

	                   for mysample in mysamples:
	                       for c in alldata[mysample][mygene]:
	                           each.append(c)

	                   print>>self.output1,mygene+'\t'+'\t'.join(each)


	        ## myoutput.close()
	        self.input.close()
	        self.output1.close()
	        self.output2.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	generate().generate()



if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
