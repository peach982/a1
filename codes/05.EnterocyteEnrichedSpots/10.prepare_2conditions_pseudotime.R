
#source activate /Lustre01/tangqianzi/anaconda3forSeurat4/envs/R4.0.2
#export HDF5_USE_FILE_LOCKING='FALSE'

library(gtools)
library(Seurat)
library(dplyr)
library(ggplot2)

input_nu<-'water'
input_depth<-5
input_locus1<-'Jejunum'
input_locus2<-'ACaecum'
#input_locus2<-'PColon'

setwd('/data/tangyuanling/STsnRNA/forQinglin20251014/monocle23_time_course/nutrient_MEs/')
data<-read.table(file=paste(input_nu,'_depth_',input_depth,'_',input_locus1,'_',input_locus2,'.txt',sep=''),sep='\t',header=T)
DMEs<-as.data.frame(data)
cp <- "orange"; names(cp) <- paste(input_nu,'_depth_',input_depth,'_',input_locus1,'_',input_locus2,sep='')
p <- DMEs %>%
	ggplot(aes(x = group, y = avg_log2FC, fill=module, color=module)) +
	geom_hline(yintercept=0, linetype='dashed', color='darkgrey') +
	geom_segment(aes(y=0, yend=avg_log2FC, x=group, xend=group), size=0.5, color='grey') +
	geom_point(size=3) +
	geom_point(shape=DMEs$shape, color='black', fill=NA, size=3) +
	scale_y_continuous(limits=c(-4, 5)) +
	scale_color_manual(values=cp, guide='none') +
	ylab(bquote("Average log"[2]~"(Fold Change)")) +
	theme(
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank()
	)
pdf(paste('diff_2conditions/',input_nu,'_depth_',input_depth,'_',input_locus1,'_',input_locus2,'.pdf',sep=''),width=6,height=3)
p+facet_wrap(~module,ncol=1)
dev.off()
