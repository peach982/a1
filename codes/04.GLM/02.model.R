##source activate /Lustre01/tangqianzi/anaconda3forSeurat4/envs/R4.0.2

library(MASS)
library(Matrix)
library(tidyverse)
library(multcomp)

setwd('/Lustre03/data/tangqianzi/Gut_unroll/outputs/ANOVAlike_analyses/')
data1<-read.table(file='01.total_counts.xls',sep='\t',header=T,row.names=1)
data2<-read.table(file='00.merged_table.xls',sep='\t',header=T,row.names=1)
data3<-read.table(file='00.info.xls',sep='\t',header=T)

iter_num<-nrow(data2)
##iter_num<-10

model_ps<-c()
depth_ps<-c()
age_ps<-c()
depth_coeffs<-c()
age_coeffs<-c()

J_A_ps<-c()
P_A_ps<-c()
P_J_ps<-c()

J_A_coeffs<-c()
P_A_coeffs<-c()
P_J_coeffs<-c()

genenames<-c()

flag<-F

for(i in 1:iter_num){

m<-i
genename<-rownames(data2)[m]

mynewdata<-cbind(data3,t(data2[m,]),data1[,1])
colnames(mynewdata)<-c(colnames(data3),'count','total')
mytotal<-as.numeric(mynewdata$total)
locus <- NULL
mynewdata$locus <- factor(mynewdata$locus)

tryCatch({
m1 <-glm.nb(count~depth+age_number+locus+log(total),data=mynewdata)},
warning = function(w){

genenames<<-append(genenames,genename)
model_ps<<-append(model_ps,'warn')
depth_ps<<-append(depth_ps,'warn')
age_ps<<-append(age_ps,'warn')
depth_coeffs<<-append(depth_coeffs,'warn')
age_coeffs<<-append(age_coeffs,'warn')
J_A_ps<<-append(J_A_ps,'warn')
P_A_ps<<-append(P_A_ps,'warn')
P_J_ps<<-append(P_J_ps,'warn')
J_A_coeffs<<-append(J_A_coeffs,'warn')
P_A_coeffs<<-append(P_A_coeffs,'warn')
P_J_coeffs<<-append(P_J_coeffs,'warn')
flag<<-T
},error = function(e){
    if(flag == F){
    genenames<<-append(genenames,genename)
    model_ps<<-append(model_ps,'error')
    depth_ps<<-append(depth_ps,'error')
    age_ps<<-append(age_ps,'error')
    depth_coeffs<<-append(depth_coeffs,'error')
    age_coeffs<<-append(age_coeffs,'error')
    J_A_ps<<-append(J_A_ps,'error')
    P_A_ps<<-append(P_A_ps,'error')
    P_J_ps<<-append(P_J_ps,'error')
    J_A_coeffs<<-append(J_A_coeffs,'error')
    P_A_coeffs<<-append(P_A_coeffs,'error')
    P_J_coeffs<<-append(P_J_coeffs,'error')
    flag<<-T}
  },finally = {
    if(flag==F){
    depth_p<-summary(m1)[[12]][2,4]
    age_p<-summary(m1)[[12]][3,4]
    depth_coeff<-summary(m1)[[12]][2,1]
    age_coeff<-summary(m1)[[12]][3,1]

    mynewdata$locus<-as.factor(mynewdata$locus)

    mult_comp <- summary(glht(m1, mcp(locus="Tukey")))

    J_A_p<-mult_comp$test$pvalues[1]  #Jejunum - ACaecum
    P_A_p<-mult_comp$test$pvalues[2]  #PColon - ACaecum
    P_J_p<-mult_comp$test$pvalues[3]  #PColon - Jejunum

    J_A_coeff<-mult_comp$test$coefficients[1]
    P_A_coeff<-mult_comp$test$coefficients[2]
    P_J_coeff<-mult_comp$test$coefficients[3]

    #m2 <-glm.nb(count~depth+age_number+locus+log(total)+locus:depth,data=mynewdata)
    #drop_in_dev <- anova(m1, m2, test = "Chisq")

    model_p<-1-pchisq(m1$deviance, m1$df.residual)

    genenames<<-append(genenames,genename)
    model_ps<<-append(model_ps,model_p)
    depth_ps<<-append(depth_ps,depth_p)
    age_ps<<-append(age_ps,age_p)
    depth_coeffs<<-append(depth_coeffs,depth_coeff)
    age_coeffs<<-append(age_coeffs,age_coeff)
    J_A_ps<<-append(J_A_ps,J_A_p)
    P_A_ps<<-append(P_A_ps,P_A_p)
    P_J_ps<<-append(P_J_ps,P_J_p)
    J_A_coeffs<<-append(J_A_coeffs,J_A_coeff)
    P_A_coeffs<<-append(P_A_coeffs,P_A_coeff)
    P_J_coeffs<<-append(P_J_coeffs,P_J_coeff)
    flag<<-T}
  })
  flag<<-F
}

model_fdrs<-model_ps
model_fdrs[model_ps!='error'&model_ps!='warn']<-p.adjust(model_ps[model_ps!='error'&model_ps!='warn'])
depth_fdrs<-depth_ps
depth_fdrs[depth_ps!='error'&depth_ps!='warn']<-p.adjust(depth_ps[depth_ps!='error'&depth_ps!='warn'])
age_fdrs<-age_ps
age_fdrs[age_ps!='error'&age_ps!='warn']<-p.adjust(age_ps[age_ps!='error'&age_ps!='warn'])
J_A_fdrs<-J_A_ps
J_A_fdrs[J_A_ps!='error'&J_A_ps!='warn']<-p.adjust(J_A_ps[J_A_ps!='error'&J_A_ps!='warn'])
P_A_fdrs<-P_A_ps
P_A_fdrs[P_A_ps!='error'&P_A_ps!='warn']<-p.adjust(P_A_ps[P_A_ps!='error'&P_A_ps!='warn'])
P_J_fdrs<-P_J_ps
P_J_fdrs[P_J_ps!='error'&P_J_ps!='warn']<-p.adjust(P_J_ps[P_J_ps!='error'&P_J_ps!='warn'])

finalresult<-cbind(model_ps,model_fdrs,depth_ps,depth_fdrs,depth_coeffs,age_ps,age_fdrs,age_coeffs,J_A_ps,J_A_fdrs,J_A_coeffs,P_A_ps,P_A_fdrs,P_A_coeffs,P_J_ps,P_J_fdrs,P_J_coeffs)
rownames(finalresult)<-genenames

colnames(finalresult)<-cbind('model_ps','model_fdrs','depth_ps','depth_fdrs','depth_coeffs','age_ps','age_fdrs','age_coeffs','J_A_ps','J_A_fdrs','J_A_coeffs','P_A_ps','P_A_fdrs','P_A_coeffs','P_J_ps','P_J_fdrs','P_J_coeffs')

write.table(finalresult,file='finalresult.xls',sep='\t',quote=F,row.names=T,col.names=T)
