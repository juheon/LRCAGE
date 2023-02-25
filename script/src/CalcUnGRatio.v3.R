
args=commandArgs(trailingOnly = T)
#args=c("/scratch/jmaeng/longreadCAGE/h1299/analysis_hg38_noalt.v2/paraclu_test/LRCAGE.BR1_2.3M.bam.CTSS.v3.consclus.bed", 
#       "/scratch/jmaeng/longreadCAGE/h1299/analysis_hg38_noalt.v2/paraclu_test/LRCAGE.BR1_2.3M.bam.CTSS.v3.consclus.bed.read_in_peak", 
#       "/scratch/jmaeng/longreadCAGE/h1299/analysis_hg38_noalt.v2/paraclu_test/LRCAGE.BR1_2.3M.bam.CTSS.v3.consclus.bed.unGread_in_peak")



#args=c("/scratch/jmaeng/longreadCAGE/h1299/pacbio/downsample.v2/3.2M_peak/LRCAGE_realn.RQ90.BR1_2.bam.consclus_3.bed",
#       "/scratch/jmaeng/longreadCAGE/h1299/pacbio/downsample.v2/3.2M_peak/LRCAGE_realn.RQ90.BR1_2.bam.consclus_3.bed.read_in_peak",
#       "/scratch/jmaeng/longreadCAGE/h1299/pacbio/downsample.v2/3.2M_peak/LRCAGE_realn.RQ90.BR1_2.bam.consclus_3.bed.unGread_in_peak",
#       "0.2",
#       "2")

ipeak=read.table(args[1], stringsAsFactors = F)
ictss=read.table(args[2], stringsAsFactors = F)
iugctss=read.table(args[3], stringsAsFactors = F)
oprefix=args[4]
gcapcutoff=as.numeric(args[5]) #gcap cutoff
gcaptag=as.numeric(args[6])    #min gcap tag count

#Count read per peak
library(dplyr)
ictss.count=aggregate( V11 ~ V4, ictss[ictss$V11>0, ], FUN=sum)
iugctss.count=aggregate( V11 ~ V4, iugctss[iugctss$V11>0, ], FUN=sum)
colnames(ictss.count)=c("V4", "readcnt")
colnames(iugctss.count)=c("V4", "unGreadcnt")

dfm=merge( ipeak, ictss.count, by="V4", all.x=T)
dfm=merge( dfm, iugctss.count, by="V4", all.x=T)
dfm[is.na(dfm$unGreadcnt), ]$unGreadcnt=0
dfm$unGratio=dfm$unGreadcnt/dfm$readcnt

vecChrom=paste0("chr", c(seq(1, 22), "X", "Y"))
dfm$V1=factor(dfm$V1, levels=vecChrom)
dfm.sort=dfm[order(dfm$V1, dfm$V2), c("V1", "V2", "V3", "V4", "V5", "V6", "readcnt", "unGreadcnt", "unGratio")]
write.table( dfm.sort, paste0( oprefix, ".unG.bed"), quote=F, sep="\t", row.names=F, col.names=F)

suffix=paste0(".unG_",gcapcutoff )
write.table( dfm.sort[ dfm.sort$unGratio>=gcapcutoff & dfm.sort$unGreadcnt>=1, ], paste0( paste0(oprefix, suffix), ".bed"), quote=F, sep="\t", row.names=F, col.names=F)


