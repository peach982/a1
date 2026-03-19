
from optparse import OptionParser
import os

class autoremove_bam:

	def __init__ (self):

	        path='/Lustre03/data/tangqianzi/Gut_unroll/outputs/ANOVAlike_analyses/'

	        #self.input1 = open('/Lustre01/tangqianzi/Gut_project/new/rawcounts/counts_'+sample+'.csv','r')
	        ##self.input1 = open('/Lustre03/data/tangqianzi/Gut_myproject_new/'+sample+'_rawcounts.txt','r')
	        #self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/pos2spotID.absolutedistance_finalgroupsforpseudobulk.xls','r')
	        ##self.input2 = open('/Lustre03/data/tangqianzi/Gut_unroll/outputs/'+sample+'/spot2group_relative_distance.'+cutoff+'.xls','r')
	        #self.input = open(path+'/info.txt','r')
	        self.input = open(path+'/finalresult.xls','r')
	        self.output = open(path+'/finalresult.genes.xls','w')
	        #self.output = open(path+'/ANOVAlike_analyses/01.total_counts.xls','w')
	        self.path = path
	        self.fdrcutoff = 0.05
	        self.folder = 'gene_lists'

	def autoremove_bam (self):

	        os.system('mkdir -p '+self.path+'/'+self.folder+'/')

                depth_pos={}
                depth_neg={}
                age_pos={}
                age_neg={}

                J_A_pos={}
                J_A_neg={}
                P_A_pos={}
                P_A_neg={}
                P_J_pos={}
                P_J_neg={}

                index2name={}
                i=0
	        for line in self.input:
	            line=line.rstrip()
	            parts=line.split('\t')
	            if i==0:
	                for j in range(0,len(parts)):
	                    index2name[str(j+1)]=parts[j]

	            else:
	                sig_keep=0
	                for j in range(1,len(parts)):
	                    if index2name[str(j)]=='model_ps':
	                        ##print parts[j]
	                        if parts[j]!='error' and parts[j]!='warn' and float(parts[j])>0.05:
	                            sig_keep=1

	                    if index2name[str(j)]=='depth_fdrs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               depth_fdrs=float(parts[j])

	                    if index2name[str(j)]=='depth_coeffs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               depth_ceoffs=float(parts[j])


	                    if index2name[str(j)]=='age_fdrs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               age_fdrs=float(parts[j])

	                    if index2name[str(j)]=='age_coeffs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               age_ceoffs=float(parts[j])



	                    if index2name[str(j)]=='J_A_fdrs':
	                        if parts[j]!='error' and parts[j]!='warn':
   	                                J_A_fdrs=float(parts[j])

	                    if index2name[str(j)]=='J_A_coeffs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                                J_A_ceoffs=float(parts[j])



	                    if index2name[str(j)]=='P_A_fdrs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               P_A_fdrs=float(parts[j])

	                    if index2name[str(j)]=='P_A_coeffs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               P_A_ceoffs=float(parts[j])

	                    if index2name[str(j)]=='P_J_fdrs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               P_J_fdrs=float(parts[j])

	                    if index2name[str(j)]=='P_J_coeffs':
	                        if parts[j]!='error' and parts[j]!='warn':
	                               P_J_ceoffs=float(parts[j])


                        if sig_keep==1:
                            if depth_fdrs<float(self.fdrcutoff) and depth_ceoffs>0:
                                    depth_pos[parts[0]]=0
                            if depth_fdrs<float(self.fdrcutoff) and depth_ceoffs<0:
                                    depth_neg[parts[0]]=0

                            if age_fdrs<float(self.fdrcutoff) and age_ceoffs>0:
                                    age_pos[parts[0]]=0
                            if age_fdrs<float(self.fdrcutoff) and age_ceoffs<0:
                                    age_neg[parts[0]]=0

                            if J_A_fdrs<float(self.fdrcutoff) and J_A_ceoffs>0:
                                    J_A_pos[parts[0]]=0
                            if J_A_fdrs<float(self.fdrcutoff) and J_A_ceoffs<0:
                                    J_A_neg[parts[0]]=0

                            if P_A_fdrs<float(self.fdrcutoff) and P_A_ceoffs>0:
                                    P_A_pos[parts[0]]=0
                            if P_A_fdrs<float(self.fdrcutoff) and P_A_ceoffs<0:
                                    P_A_neg[parts[0]]=0

                            if P_J_fdrs<float(self.fdrcutoff) and P_J_ceoffs>0:
                                    P_J_pos[parts[0]]=0
                            if P_J_fdrs<float(self.fdrcutoff) and P_J_ceoffs<0:
                                    P_J_neg[parts[0]]=0

	            i+=1

                print>>self.output,'gene_number'+'\t'+'genenames'
                print>>self.output,'depth_pos'+'\t'+str(len(depth_pos))+'\t'+','.join(depth_pos)
                print>>self.output,'depth_neg'+'\t'+str(len(depth_neg))+'\t'+','.join(depth_neg)
                print>>self.output,'age_pos'+'\t'+str(len(age_pos))+'\t'+','.join(age_pos)
                print>>self.output,'age_neg'+'\t'+str(len(age_neg))+'\t'+','.join(age_neg)

                print>>self.output,'J_A_pos'+'\t'+str(len(J_A_pos))+'\t'+','.join(J_A_pos)
                print>>self.output,'J_A_neg'+'\t'+str(len(J_A_neg))+'\t'+','.join(J_A_neg)
                print>>self.output,'P_A_pos'+'\t'+str(len(P_A_pos))+'\t'+','.join(P_A_pos)
                print>>self.output,'P_A_neg'+'\t'+str(len(P_A_neg))+'\t'+','.join(P_A_neg)
                print>>self.output,'P_J_pos'+'\t'+str(len(P_J_pos))+'\t'+','.join(P_J_pos)
                print>>self.output,'P_J_neg'+'\t'+str(len(P_J_neg))+'\t'+','.join(P_J_neg)

                handle=open(self.path+'/'+self.folder+'/depth_pos.txt','w')
                for c in depth_pos:
                    print>>handle,c
                handle.close()

                handle=open(self.path+'/'+self.folder+'/depth_neg.txt','w')
                for c in depth_neg:
                    print>>handle,c
                handle.close()

                handle=open(self.path+'/'+self.folder+'/age_pos.txt','w')
                for c in age_pos:
                    print>>handle,c
                handle.close()

                handle=open(self.path+'/'+self.folder+'/age_neg.txt','w')
                for c in age_neg:
                    print>>handle,c
                handle.close()


                handle=open(self.path+'/'+self.folder+'/J_A_pos.txt','w')
                for c in J_A_pos:
                    print>>handle,c
                handle.close()

                handle=open(self.path+'/'+self.folder+'/J_A_neg.txt','w')
                for c in J_A_neg:
                    print>>handle,c
                handle.close()

                handle=open(self.path+'/'+self.folder+'/P_A_pos.txt','w')
                for c in P_A_pos:
                    print>>handle,c
                handle.close()

                handle=open(self.path+'/'+self.folder+'/P_A_neg.txt','w')
                for c in P_A_neg:
                    print>>handle,c
                handle.close()


                handle=open(self.path+'/'+self.folder+'/P_J_pos.txt','w')
                for c in P_J_pos:
                    print>>handle,c
                handle.close()

                handle=open(self.path+'/'+self.folder+'/P_J_neg.txt','w')
                for c in P_J_neg:
                    print>>handle,c
                handle.close()

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
