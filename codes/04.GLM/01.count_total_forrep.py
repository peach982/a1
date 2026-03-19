
from optparse import OptionParser

class autoremove_bam:

	def __init__ (self):

	        path='/Lustre03/data/tangqianzi/Gut_unroll/outputs/'

	        #self.input1 = open('/Lustre01/tangqianzi/Gut_project/new/rawcounts/counts_'+sample+'.csv','r')
	        ##self.input1 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/'+sample+'_rawcounts.txt','r')
	        #self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.absolutedistance_finalgroupsforpseudobulk.xls','r')
	        ##self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/spot2group_relative_distance.'+cutoff+'.xls','r')
	        #self.input = open(path+'/info.txt','r')
	        self.input = open(path+'/ANOVAlike_analyses/00.merged_table.xls','r')
	        self.output = open(path+'/ANOVAlike_analyses/01.total_counts.xls','w')
	        self.path = path

	def autoremove_bam (self):

	        i=0
	        allcounts={}
	        keepnames=[]
	        for line in self.input:
	            line=line.rstrip()
	            parts=line.split('\t')
	            if i==0:
	                keepnames=parts
	                for j in range(0,len(parts)):
	                    allcounts[parts[j]]=0
	            else:
	                for j in range(1,len(parts)):
	                    allcounts[keepnames[j-1]]+=int(parts[j])

                    i+=1


                print>>self.output,'total_counts'
                for myname in keepnames:
                    print>>self.output,myname+'\t'+str(allcounts[myname])

	        self.input.close()
	        self.output.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	autoremove_bam().autoremove_bam()


if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
