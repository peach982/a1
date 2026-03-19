
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time

#export PATH=/Lustre01/tangqianzi/software/anaconda2new/bin/:$PATH
#/lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
#/Lustre03/data/tangqianzi/Gut_timepoint_bulkRNAseq/ensembl-pig-add-MYH24.gtf

class generate:

	def __init__ (self,path='/Lustre03/data/tangqianzi/Gut_myproject_new/',name='2J_2',cutoff='6'):

		self.path = path+'/'
		self.input1 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+name+'/finalresult_pseudobulk_relativedistance.'+cutoff+'.xls','r')
		self.input2 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/limmaresultsfinal/pig2human.txt','r')
		self.input6 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/limmaresultsfinal/gene2types.xls','r')
		self.output = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+name+'/finalresult_pseudobulk_relativedistance.'+cutoff+'.humanids.xls','w')
		#self.name = name

	def generate (self):

	        #handleout1 = open(self.path+'/'+newname+'.'+myname1+'.forDAVID.txt','w')
	        #handleout2 = open(self.path+'/'+newname+'.'+myname2+'.forDAVID.txt','w')


	        piggenekeep={}
	        i=0
	        for line in self.input6:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i!=0 and parts[1]=='protein_coding':
	                       piggenekeep[parts[0]]=0

	               i+=1

	        i=0
	        ##mynames_1=[]
	        ##mynames_2=[]
	        #headline=''
	        mydata={}
                exclude={}
	        for line in self.input1:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i==0:
	                       ##headline=line
	                       print>>self.output,'humanname'+'\t'+'pigname'+'\t'+'\t'.join(parts[:-2])+'\t'+'type'

	               elif i!=0:
	                       if not(parts[0].startswith('ENSSSCG')):
	                               print>>self.output,parts[0]+'\t'+'\t'.join(parts[:-2])+'\t'+'direct'
	                               exclude[parts[0]]=0

	                       else:
	                               mydata[parts[0]]='\t'.join(parts[:-2])








	               i+=1



                i=0
                alldata={}
                human2pig={}
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

                            print parts[6]

                            alldata[pigname].append([int(parts[8]),float(parts[5]),float(parts[6]),float(parts[7]),parts[3],parts[4]])

                            if parts[3] not in human2pig:
                                    human2pig[parts[3]]=[]

                            human2pig[parts[3]].append(pigname)


                       i+=1

                count=1
                myconvert={}

                count_not=0

                for myname in alldata:
                       alldata[myname].sort()
                       alldata[myname].reverse()

                       #if not(re.search('ssc-',alldata[myname][0][2])) and not(re.search('MT-',alldata[myname][0][3])) and alldata[myname][0][2]!=alldata[myname][0][3] and alldata[myname][0][2]!='' and alldata[myname][0][3]!='':
                            ##print alldata[myname][0][2],alldata[myname][0][3]
                            #count+=1

                       if len(alldata[myname])>1 and alldata[myname][0][0]==alldata[myname][1][0] and alldata[myname][0][1]==alldata[myname][1][1] and alldata[myname][0][2]==alldata[myname][1][2] and alldata[myname][0][3]==alldata[myname][1][3]:
                            count_not+=1

                       if myname not in exclude:
                            myconvert[myname]=[alldata[myname][0][4],alldata[myname][0][5]]


                print count

                print "count_not",count_not

                finalkeep={}
                for pigname in mydata:
                        if pigname in alldata:
                                myline=mydata[pigname]
                                humanname=myconvert[pigname][0]

                                if humanname not in finalkeep:
                                        finalkeep[humanname]=[]

                                finalkeep[humanname].append([myline,myconvert[pigname][1]])

                for humanname in finalkeep:
                        if humanname!='' and (humanname not in exclude):
                                if len(finalkeep[humanname])==1:

                                        keeptypes=finalkeep[humanname][0][1]
                                        print>>self.output,humanname+'\t'+finalkeep[humanname][0][0]+'\t'+keeptypes
                                else:
                                        mynumele=len(finalkeep[humanname][0][0].split('\t'))-1
                                        myeach_sum=[0]*mynumele
                                        myeach_ave=[0]*mynumele

                                        mynumgene=len(finalkeep[humanname])

                                        keeppignames=[]
                                        keeptypes=[]

                                        for m in range(0,mynumgene):

                                                keeppignames.append(finalkeep[humanname][m][0].split('\t')[0])
                                                keeptypes.append(finalkeep[humanname][m][1])
                                                eacheach=finalkeep[humanname][m][0].split('\t')[1:]
                                                for k in range(0,len(eacheach)):
                                                        myeach_sum[k]+=float(eacheach[k])

                                        for k in range(0,mynumele):
                                                myeach_ave[k]=str(myeach_sum[k]/float(mynumgene))

                                        print>>self.output,humanname+'\t'+','.join(keeppignames)+'\t'+'\t'.join(myeach_ave)+'\t'+','.join(keeptypes)

		self.input1.close()
		self.input2.close()
		self.input6.close()
		self.output.close()

	        #handleout1.close()
	        #handleout2.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	#generate().generate()

	generate(cutoff=pathandfiles[0],name=pathandfiles[1]).generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
