
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

class autoremove_bam:

	def __init__ (self,sample='2J_1'):

		#self.pathin = pathin
		self.input = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/tissue_positions/'+sample+'/tissue_positions_list.csv','r')

                #change this every time
		self.input2 = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/2J_1-XYZ-line-0115.csv','r')

	        self.output = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/outputs/'+sample+'/orderX.R','w')
	        #self.output2 = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/outputs/'+sample+'/simulate_spots.sameasothers.XYZ.X.test.xls','w')
	        #self.output3 = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/outputs/'+sample+'/pos2spotID.test.xls','w')
	        self.output4 = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/outputs/'+sample+'/orderX.xls','w')
	        self.output5 = open('/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/outputs/'+sample+'/orderX.exclude.xls','w')

                self.sample = sample

	def autoremove_bam (self):

                def line_intersection(line1, line2):
                        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
                        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

                        def det(a, b):
                                return a[0] * b[1] - a[1] * b[0]

                        div = det(xdiff, ydiff)
                        if div == 0:
                                raise Exception('lines do not intersect')

                        d = (det(*line1), det(*line2))
                        x = det(d, xdiff) / div
                        y = det(d, ydiff) / div
                        return x, y

	        inside_spots=[]
	        outside_spots=[]
	        in_spots_old=[]
	        in_spots=[]

	        i=0
	        for line in self.input2:
	               line=line.rstrip()
	               parts=line.split(',')
	               if i!=0 and len(parts)>=2:
	                       if parts[1]=='XX':
	                               inside_spots.append(parts[0])
	                       elif parts[1]=='X_Y':
	                               outside_spots.append(parts[0])
	                       elif parts[1]=='X_inside':
	                       #elif parts[1]=='X':
	                               in_spots_old.append(parts[0])



	               i+=1

	        print len(in_spots)

	        print>>self.output,'setwd("/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/outputs/'+self.sample+'/")'
	        print>>self.output,'pdf("orderX.pdf",h=12,w=10,bg="transparent")'
	        xlim2=math.sqrt(3)*77+1
	        ## xlim2=127+1

	        print>>self.output,'layout(matrix(1:2,nrow=2),c(10),c(5,1))'
	        print>>self.output,'par(mar=c(2,2,2,2))'
	        print>>self.output,'plot(c(1,1),type="n",xlab="",ylab="",xlim=c(-1,'+str(xlim2)+'),ylim=c(-1,'+str(xlim2)+'),axes=F)'
	        print>>self.output,'library(plotrix)'

                h=0
                id2pos={}
                pos2id={}
                pos2id_2={}

                allallallx=[]
	        for line in self.input:
	               line=line.rstrip()
	               parts=line.split(',')

	               #if h!=0:
	               if 1:
	                       myx=int(parts[3])
	                       myy=xlim2-int(parts[2])*math.sqrt(3)

	                       allallallx.append(myx)




	                       mysize=(55/50.0)/2.0

	                       print>>self.output,'draw.circle('+str(myx)+','+str(myy)+','+str(mysize)+',border="darkgreen",lwd=0.5)'
	                       ## DrawCircle(x = 0, y = x, r.out = 1

	                       #print>>self.output,'draw.circle('+str(myx2)+','+str(myy2)+','+str(mysize)+',border="darkred",lwd=0.5)'

	                       id2pos[parts[0]]=[int(parts[2]),int(parts[3])]

	                       pos2id[str(myx)+'_'+str(myy)]=parts[0]
	                       pos2id_2[parts[2]+'_'+parts[3]]=parts[0]



	               h+=1



	        #exclude_outside={}
	        #exclude_inside={}

                allrestspots_out_old=[]
                allkeepspots_out_old=[]

                allrestspots_in=[]
                allkeepspots_in=[]


                print len(outside_spots)

                allallallx.sort()
                allallallx.reverse()

                oldc=127
                for c in allallallx:
                        if c!=oldc:
                            #print c
                            pass



                county0=0
                countx0=0
                countxmax=0

                alltestx2=[]
	        for myname2 in outside_spots:
	               y2=xlim2-id2pos[myname2][0]*math.sqrt(3)
	               x2=id2pos[myname2][1]

	               alltestx2.append(x2)

	               if x2==max(allallallx):
	                       #startx_out=x2
	                       #starty_out=y2
	                       allkeepspots_out_old.append(myname2)
	                       countxmax+=1

                       else:
                               #not sure correct or not???
                               #startx_out=0
	                       allrestspots_out_old.append(myname2)


	        print max(alltestx2),"max(alltestx2)"

	        print countx0,"countx0"
	        print county0,"county0"
	        print countxmax,"countxmax"




	        if countxmax>=2:

	               eacheach=[]
	               for myname2 in allkeepspots_out_old:
	                       y2=xlim2-id2pos[myname2][0]*math.sqrt(3)
	                       x2=id2pos[myname2][1]

	                       eacheach.append([y2,myname2])
	               eacheach.sort()


	               mykeep1name=eacheach[0][1]
	               myremove1name=[]

	               for m in range(1,len(eacheach)):
	                       myremove1name.append(eacheach[m][1])



	               allkeepspots_out=[]
	               allrestspots_out=[]

	               allkeepspots_out.append(mykeep1name)

	               for c in allrestspots_out_old:
	                       allrestspots_out.append(c)

	               for c in myremove1name:
	                       allrestspots_out.append(c)

	               startx_out=id2pos[mykeep1name][1]
	               starty_out=xlim2-id2pos[mykeep1name][0]*math.sqrt(3)




	        for myname2 in inside_spots:
	               y2=xlim2-id2pos[myname2][0]*math.sqrt(3)
	               x2=id2pos[myname2][1]

	               if y2==xlim2:
	                       startx_in=0
	                       starty_in=y2
	                       allkeepspots_in.append(myname2)

                       else:
	                       allrestspots_in.append(myname2)


                def get_order(start_x,start_y,allrestspots,keepspots):

                       if len(allrestspots)>0:

                            alldists=[]
	                    for i in range(0,len(allrestspots)):
	                           myy=xlim2-id2pos[allrestspots[i]][0]*math.sqrt(3)
	                           myx=id2pos[allrestspots[i]][1]
	                           mydist=math.sqrt((start_x-myx)**2+(start_y-myy)**2)

	                           if mydist>0:

	                                   alldists.append([mydist,allrestspots[i]])

	                    alldists.sort()

	                    allrestspots.remove(alldists[0][1])
	                    keepspots.append(alldists[0][1])

	                    start_y=xlim2-id2pos[alldists[0][1]][0]*math.sqrt(3)
	                    start_x=id2pos[alldists[0][1]][1]

	                    get_order(start_x,start_y,allrestspots,keepspots)

	               else:

	                    return 1

                #outside_spots_ordered=get_order(startx_out,starty_out,allrestspots_out,allkeepspots_out)
                correctness=get_order(startx_out,starty_out,allrestspots_out,allkeepspots_out)

                ##correctness2=get_order(startx_in,starty_in,allrestspots_in,allkeepspots_in)

                print len(allkeepspots_in),"test0"

                #print len(outside_spots_ordered),"test1"
                print len(inside_spots),"test2"




	        for m in outside_spots:
	               y1=xlim2-id2pos[m][0]*math.sqrt(3)
	               x1=id2pos[m][1]

                       print>>self.output,'draw.circle('+str(x1)+','+str(y1)+','+str(mysize)+',border="purple",col="purple",lwd=0.5)'

	        for m in inside_spots:
	               y1=xlim2-id2pos[m][0]*math.sqrt(3)
	               x1=id2pos[m][1]

                       print>>self.output,'draw.circle('+str(x1)+','+str(y1)+','+str(mysize)+',border="yellow",col="yellow",lwd=0.5)'

                       if m=='AACAGGATGGGCCGCG-1':
	                       print "test0"


                if 1:
                       mycount=1
                       for i in range(0,len(allkeepspots_out)):
                            startname=allkeepspots_out[i]
                            #endname=allkeepspots_out[i+1]
                            startp_y=xlim2-id2pos[startname][0]*math.sqrt(3)
                            #endp_x=id2pos[endname][0]*math.sqrt(3)

                            startp_x=id2pos[startname][1]

                            print>>self.output,'draw.circle('+str(startp_x)+','+str(startp_y)+','+str(mysize)+',border="purple",lwd=0.5)'
                            print>>self.output,'text('+str(startp_x)+','+str(startp_y)+','+str(mycount)+',cex=0.3)'

                            if mycount==67:
                                    print startname
                                    myexclude_name=startname
                                    print>>self.output5,startname

                            ##print>>self.output4,startname

                            mycount+=1

                #hasdata=[]

                #for c in exclude_inside:
                       #myx1=id2pos[c][0]*math.sqrt(3)
                       #myy1=xlim2-id2pos[c][1]
                       #print>>self.output,'draw.circle('+str(myx1)+','+str(myy1)+','+str(mysize)+',col="orange",lwd=0.5)'

                #for c in exclude_outside:
                       #myx1=id2pos[c][0]*math.sqrt(3)
                       #myy1=xlim2-id2pos[c][1]
                       #print>>self.output,'draw.circle('+str(myx1)+','+str(myy1)+','+str(mysize)+',col="blue",lwd=0.5)'


                #for c in in_spots_old:
                       #if (c not in exclude_inside) and (c not in exclude_outside):
                               #in_spots.append(c)


                for c in in_spots_old:
                       in_spots.append(c)

                #print len(exclude_inside)
                #print len(exclude_outside)
                print len(in_spots_old)
                print len(in_spots)

                ##print>>self.output3,'start_id'+'\t'+'end_id'+'\t'+'distance'

                #print>>self.output2,'start_id'+'\t'+'end_id'+'\t'+'distance'

                print>>self.output4,'purple_out'+'\t'+'yellow_inside'

                countall=0
                count_wrong=0
                count_print=0
	        for p in outside_spots:
	               #if m in inside_spots:
	                       #print "test",m

	               #if p=='AACAGGATGGGCCGCG-1':
	                       #print "test"


	               #if countall!=28:
	                       #countall+=1
	                       #continue

	               if p==myexclude_name:
	                       continue

	               #print countall
	               y1=xlim2-id2pos[p][0]*math.sqrt(3)
	               x1=id2pos[p][1]

	               alldists=[]
	               for n in inside_spots:
	                       #if n!=m:
	                       if 1:

	                               y2=xlim2-id2pos[n][0]*math.sqrt(3)
	                               x2=id2pos[n][1]

	                               mydist=math.sqrt((x1-x2)**2+(y1-y2)**2)

	                               alldists.append([mydist,x2,y2,n])

	               alldists.sort()

	               #from small to large

                       sighas=0
                       count1=0
	               for parts in alldists:

	                       #print count1
	                       sigcross=0
	                       sig_equal=0
	                       sigtest=0
	                       #sig_equaltest=0

	                       myfinaly=parts[2]
	                       myfinalx=parts[1]

	                       #outside spot
	                       myfinalspot=parts[3]

                               all_higher=[]
                               all_lower=[]

	                       myx1=x1
	                       myx2=myfinalx

	                       myy1=y1
	                       myy2=myfinaly

                               #sig_equal3=0
                               #inner line

                               mytype='NA'
	                       for k in outside_spots:

	                               #if str(k)!=str(p):
	                               if 1:
	                                     y3=xlim2-id2pos[k][0]*math.sqrt(3)
	                                     x3=id2pos[k][1]

	                                     if x3==myx1 and y3==myy1:
	                                           continue

                                             if float(myfinalx-myx1)!=0 and float(myfinaly-myy1)!=0:
	                                           mval=(myfinaly-myy1)/float(myfinalx-myx1)

	                                           y3new=(x3-myx1)*mval+myy1


	                                           #distance from point to line

	                                           B=(myfinalx-myx1)
	                                           A=-(myfinaly-myy1)
	                                           C=myx1*myfinaly-myfinalx*myy1

	                                           fenzi=abs(A*x3+B*y3+C)
	                                           fenmu=math.sqrt(A**2+B**2)

	                                           mydistance2=fenzi/float(fenmu)

                                                   #print mydistance2,"hi"
                                                   #if mydistance2==0:
	                                           #if mydistance2<2.0*55/100.0:
	                                           if mydistance2<=0.2:
	                                           #if mydistance2==0:
	                                                   #print mydistance2
	                                                   #print>>self.output,'lines(c('+str(myx3)+','+str(myx4)+'),c('+str(myy3)+','+str(myy4)+'),col="purple",lwd=0.5)'
	                                                   #print>>self.output,'draw.circle('+str(x3)+','+str(y3)+','+str(0.1)+',border="purple",col="purple",lwd=0.5)'
	                                                   if y3>=min(myy1,myfinaly) and y3<=max(myy1,myfinaly) and x3>=min(myx1,myfinalx) and x3<=max(myx1,myfinalx):
	                                                           #print>>self.output,'draw.circle('+str(x3)+','+str(y3)+','+str(0.3)+',border="orange",col="orange",lwd=0.5)'
	                                                           sig_equal=1
	                                                           #print>>self.output,'lines(c('+str(myx1)+','+str(myx2)+'),c('+str(myy1)+','+str(myy2)+'),col="orange",lwd=0.5)'
	                                                           ##print>>self.output,'draw.circle('+str(x3)+','+str(y3)+','+str(mysize)+',border="orange",col="orange",lwd=0.5)'

	                                                           #sigtest=1
	                                                           break
	                                                   else:
	                                                           #pass
	                                                           if y3new<y3:
	                                                                   all_higher.append([x3,y3])
	                                                           elif y3new>y3:
	                                                                   all_lower.append([x3,y3])



                                                   else:
	                                                   if y3new<y3:
	                                                           all_higher.append([x3,y3])
	                                                   elif y3new>y3:
	                                                           all_lower.append([x3,y3])

	                                     elif float(myfinalx-myx1)!=0 and float(myfinaly-myy1)==0:
	                                           if y3>myy1:
	                                                   all_higher.append([x3,y3])
	                                           elif y3<myy1:
	                                                   all_lower.append([x3,y3])
	                                           elif y3==myy1:
                                                           if (x3>=min(myx1,myfinalx) and x3<=max(myx1,myfinalx)):
	                                                           sig_equal=1
	                                                           break

	                                     elif float(myfinalx-myx1)==0 and float(myfinaly-myy1)!=0:
	                                           #y on this line can be anyvalue
	                                           #x on this line has fixed value x1, myfinalx
	                                           if x3>myx1:
	                                                   all_higher.append([x3,y3])
	                                           elif x3<myx1:
	                                                   all_lower.append([x3,y3])
	                                           elif x3==myx1:
	                                                   if (y3>=min(myy1,myfinaly) and y3<=max(myy1,myfinaly)):
	                                                           sig_equal=1
	                                                           break


	                       #if (sig_higher==1 and sig_lower==1) or (sig_equal==1):
                               #print len(all_higher)
                               #print len(all_lower)


                               all_higher.sort()
                               all_lower.sort()

                               #print "sigtest",sigtest
                               #print "siq_equal",sig_equal

                               #if sigtest!=sig_equal:
                                       #print "warning",sigtest,sig_equal
                                       #pass

                               if 1:
                                       if sig_equal!=0:
                                            sigcross=1
                                            #if str(p)=='TTGTTCTAGATACGCT-1':
                                                    #mytype='equal'

                                       elif len(all_higher)!=0 and len(all_lower)!=0:

	                                    for c in all_higher:
	                                           myx3=c[0]
	                                           myy3=c[1]


	                                           for d in all_lower:

	                                                   myx4=d[0]
	                                                   myy4=d[1]

	                                                   #if abs(myx4-myx3)>1*math.sqrt(3):

	                                                   #if math.sqrt((myx4-myx3)**2+(myy4-myy3)**2)>3.4641016151377544:
	                                                   if math.sqrt((myx4-myx3)**2+(myy4-myy3)**2)>(2+0.2):
	                                                           continue

	                                                   #xcross= -((myx1*myy2-myx2*myy1)*(myx4-myx3)-(myx3*myy4-myx4*myy3)*(myx2-myx1))/((myy1-myy2)*(myx4-myx3)-(myy3-myy4)*(myx2-myx1))


                                                           if float(myfinalx-x1)!=0 and float(myfinaly-y1)!=0:

	                                                           if (myx4-myx3)!=0:

	                                                                   #m1=myx3*(myy4-myy3)/(myx4-myx3)-myx1*(myy2-myy1)/(myx2-myx1)
	                                                                   #m2=(myy4-myy3)/(myx4-myx3)-(myy2-myy1)/(myx2-myx1)

	                                                                   #xcross=m1/m2
	                                                                   #ycross=(xcross-myx1)*(myy2-myy1)/(myx2-myx1)+myy1

	                                                                   mycoeff1=float(myy2-myy1)/float(myx2-myx1)
	                                                                   mycoeff2=float(myy4-myy3)/float(myx4-myx3)

	                                                                   #print mycoeff1,mycoeff2

	                                                                   if mycoeff1!=mycoeff2 and myy3!=myy4:

	                                                                           xcross = line_intersection([[myx1,myy1],[myx2,myy2]],[[myx3,myy3],[myx4,myy4]])[0]
	                                                                           ycross = line_intersection([[myx1,myy1],[myx2,myy2]],[[myx3,myy3],[myx4,myy4]])[1]

	                                                                           if xcross>=min(x1,myfinalx) and xcross<=max(x1,myfinalx) and ycross>=min(y1,myfinaly) and ycross<=max(y1,myfinaly) and xcross>=min(myx3,myx4) and xcross<=max(myx3,myx4) and ycross>=min(myy3,myy4) and ycross<=max(myy3,myy4):
	                                                                                   sigcross=1

                                                                                           #print>>self.output,'lines(c('+str(myx1)+','+str(myx2)+'),c('+str(myy1)+','+str(myy2)+'),col="blue",lwd=0.5)'
                                                                                           #print>>self.output,'lines(c('+str(myx3)+','+str(myx4)+'),c('+str(myy3)+','+str(myy4)+'),col="blue",lwd=0.5)'






	                                                                                   mytype='type0'
	                                                                                   mysize2=(55/50.0)/6.0

	                                                                   elif mycoeff1!=mycoeff2 and myy3==myy4:

	                                                                           ycross=myy3
	                                                                           xcross=(ycross-myy1)*(myx2-myx1)/(myy2-myy1)+myx1

	                                                                           if xcross>=min(x1,myfinalx) and xcross<=max(x1,myfinalx) and ycross>=min(y1,myfinaly) and ycross<=max(y1,myfinaly) and xcross>=min(myx3,myx4) and xcross<=max(myx3,myx4) and ycross>=min(myy3,myy4) and ycross<=max(myy3,myy4):
	                                                                                   sigcross=1
	                                                                                   mytype='type1'
	                                                                                   break

	                                                                           #if str(p)=='TTGTTCTAGATACGCT-1':
	                                                                                   #print>>self.output,'draw.circle('+str(xcross)+','+str(ycross)+','+str(mysize2)+',border="purple",col="purple",lwd=0.5)'
	                                                                                   #print>>self.output,'lines(c('+str(myx3)+','+str(myx4)+'),c('+str(myy3)+','+str(myy4)+'),col="purple",lwd=0.5)'

	                                                                           #break

	                                                           elif (myx4-myx3)==0:
	                                                                   #print "imp1"
	                                                                   #pass

                                                                           xcross=myx3
	                                                                   ycross=(xcross-myx1)*(myy2-myy1)/(myx2-myx1)+myy1


	                                                                   if xcross>=min(x1,myfinalx) and xcross<=max(x1,myfinalx) and ycross>=min(y1,myfinaly) and ycross<=max(y1,myfinaly) and xcross>=min(myx3,myx4) and xcross<=max(myx3,myx4) and ycross>=min(myy3,myy4) and ycross<=max(myy3,myy4):
	                                                                           sigcross=1
	                                                                           mytype='type1'
	                                                                           break


	                                                   elif float(myfinalx-x1)==0 and float(myfinaly-y1)!=0:
	                                                           if (myx4-myx3)!=0:

	                                                                   xcross=myx1
	                                                                   ycross= (xcross-myx3)*(myy4-myy3)/(myx4-myx3)+myy3

	                                                                   if xcross>=min(x1,myfinalx) and xcross<=max(x1,myfinalx) and ycross>=min(y1,myfinaly) and ycross<=max(y1,myfinaly) and xcross>=min(myx3,myx4) and xcross<=max(myx3,myx4) and ycross>=min(myy3,myy4) and ycross<=max(myy3,myy4):
	                                                                           sigcross=1
	                                                                           mytype='type2'
	                                                                           break

	                                                           else:
	                                                                   print "imp2"


	                                                   elif float(myfinalx-x1)!=0 and float(myfinaly-y1)==0:
	                                                           if (myx4-myx3)!=0:

	                                                                   ycross=myy1
	                                                                   xcross=(ycross-myy3)*(myx4-myx3)/(myy4-myy3)+myx3

	                                                                   if xcross>=min(x1,myfinalx) and xcross<=max(x1,myfinalx) and ycross>=min(y1,myfinaly) and ycross<=max(y1,myfinaly) and xcross>=min(myx3,myx4) and xcross<=max(myx3,myx4) and ycross>=min(myy3,myy4) and ycross<=max(myy3,myy4):
	                                                                           sigcross=1
	                                                                           mytype='type3'
	                                                                           break

	                                                           elif myx4==myx3:
	                                                                   if myx3>=min(x1,myfinalx) and myx3<=max(x1,myfinalx) and ((myy3>=y1 and myy4<=y1) or (myy3<=y1 and myy4>=y1)):
	                                                                           sigcross=1
	                                                                           mytype='type4'
	                                                                           break


	                                           if sigcross==1:
	                                                   break

	                       if sig_equal!=0 and sigcross==0:
	                               print sigcross,"mytest"

	                       if sigcross==0:
	                               #else:

	                               veryfinalx=myfinalx
	                               veryfinaly=myfinaly
	                               veryfinalspot=myfinalspot
	                               sighas=1
	                               break

	                       elif sig_equal==0 and sigcross==1:

	                               all_high_dists=[]
	                               for c in all_higher:
	                                       all_high_dists.append(math.sqrt((myx1-c[0])**2+(myy1-c[1])**2))

	                               all_low_dists=[]
	                               for c in all_lower:
	                                       all_low_dists.append(math.sqrt((myx1-c[0])**2+(myy1-c[1])**2))


                                       mydist_x1_myx3=math.sqrt((myx1-myx3)**2+(myy1-myy3)**2)
	                               mydist_x1_myx4=math.sqrt((myx1-myx4)**2+(myy1-myy4)**2)

	                               #print mydist_x1_myx3,mydist_x1_myx4


	                               #if mydist_x1_myx3==min(all_high_dists) and mydist_x1_myx4==min(all_low_dists):
	                               if mydist_x1_myx3==2 and mydist_x1_myx4==2:
	                                       print "hi"
	                                       print mydist_x1_myx3,mydist_x1_myx4
	                                       veryfinalx=myfinalx
	                                       veryfinaly=myfinaly
	                                       veryfinalspot=myfinalspot
	                                       sighas=1
	                                       break

	                       else:
	                               pass

	                               #if str(p)=='TTGTTCTAGATACGCT-1':
	                               #        print>>self.output,'lines(c('+str(x1)+','+str(myfinalx)+'),c('+str(y1)+','+str(myfinaly)+'),col="orange",lwd=0.5)'
	                               #        print>>self.output,'text(c('+str(myfinalx)+','+str(myfinaly)+'),labels="'+mytype+'",cex=0.3)'



	                       count1+=1

                       countall+=1



                       #if countall>=500:
                            #break

                       #print countall

	               #hasdata.append(alldists[0][3])

	               #print str(p),"hi0"



	               if sighas==0:
	                       #print "wrong"
	                       count_wrong+=1
	                       #print "hi",str(p)

	                       #if str(p)=='TTGTTCTAGATACGCT-1':
	                               #print>>self.output,'draw.circle('+str(myx1)+','+str(myy1)+','+str(mysize)+',col="orange",lwd=0.5)'

	                       #print>>self.output,'draw.circle('+str(x1)+','+str(y1)+','+str(mysize)+',col="pink",lwd=0.5)'


                               pass

	               else:
	                       print>>self.output,'lines(c('+str(x1)+','+str(veryfinalx)+'),c('+str(y1)+','+str(veryfinaly)+'),col="black",lwd=0.5)'

	                       endid=pos2id[str(x1)+'_'+str(y1)]
	                       startid=pos2id[str(veryfinalx)+'_'+str(veryfinaly)]


	                       mydist_final=math.sqrt((x1-veryfinalx)**2+(y1-veryfinaly)**2)

	                       #print>>self.output3,startid+'\t'+endid+'\t'+str(mydist_final)+'\t'+'X'

	                       #print>>self.output2,startid+'\t'+endid+'\t'+'X'

	                       count_print+=1

	                       print>>self.output4,endid+'\t'+startid

	               #print>>self.output,'lines(c('+str(x1)+','+str(myfinalx)+'),c('+str(y1)+','+str(myfinaly)+'),col="black",lwd=0.5)'



                print>>self.output,'box()'
                print>>self.output,'par(mar=c(2,2,2,2))'
                print>>self.output,'plot(c(1,1),type="n",xlab="",ylab="",xlim=c(-1,'+str(xlim2)+'),ylim=c(-1,1),axes=F)'

                print>>self.output,'lines(x=c(60,70),y=c(0,0))'
                print>>self.output,'text(x=c(65),y=c(-0.5),labels="500 um")'

	        print>>self.output,'dev.off()'

	        print "count_wrong",count_wrong
	        print "count_print",count_print

		self.input.close()
		self.input2.close()
		#self.input3.close()
		self.output.close()
		#self.output2.close()
		#self.output3.close()
		self.output4.close()
		self.output5.close()


		os.system('Rscript '+'/Users/hsl/Desktop/Desktop20220118/mostrecent_codes/forKong/Gut_unroll/unroll20220108/markedXYZpots/outputs/'+self.sample+'/orderX.R')



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
