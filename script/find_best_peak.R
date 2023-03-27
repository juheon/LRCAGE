args=commandArgs(trailingOnly = T)
library(dplyr)

df.readcnt=read.table(args[1], stringsAsFactors = F, header=T)
df.readcnt_pertxpeak=read.table(args[2], stringsAsFactors = F, header=T)
df.libsize=read.table(args[3], stringsAsFactors = F, header=T)

cutoff=as.numeric(args[4])

dfm=merge(df.readcnt, df.readcnt_pertxpeak, by=c("dataset", "annot_transcript_id"), all.x=T)
dfm$readcount_link=ifelse(is.na(dfm$readcount_link), 0, dfm$readcount_link)


dfm=merge( dfm, df.libsize[, c("dataset", "libsize")], by="dataset")
dfm$TPM=round( (dfm$readcount*10^6)/dfm$libsize, 2)
dfm$peak_ratio=dfm$readcount_link/dfm$readcount



#Transcript with max TPM is chosen across dataset
dfm<-dfm %>% group_by( annot_transcript_id ) %>% dplyr::mutate( IsRepresent_peaktx=ifelse( TPM==max(TPM), "O", "X" ) ) 
dfm.onepeak<-dfm[dfm$IsRepresent_peaktx=="O",!colnames(dfm) %in% c("readcount", "readcount_link", "libsize", "IsRepresent_peaktx") ] %>% group_by( annot_transcript_id ) %>% dplyr::slice_sample(n=1) 


dfms.onepeak<-dfm.onepeak %>% tidyr::spread( dataset, TPM, fill=0);
dfms.onepeak$IsConfident=ifelse(dfms.onepeak$PassMinCount=="O" & 
                                  !is.na(dfms.onepeak$peakid) & dfms.onepeak$peak_ratio>=cutoff, "O", "X")
write.table(dfms.onepeak, args[5], quote=F, sep="\t", col.names = T, row.names = F)

