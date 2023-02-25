args=commandArgs(trailingOnly=T)
library(dplyr)


df=read.table(args[1], sep="\t", stringsAsFactors = F)

#V4: peak name
#V13: dataset
#V14: transcript
#V10: read name

df.read_per_txpeak<-df %>% group_by( V1, V2, V3, V13, V14, V4 ) %>% dplyr::summarise( count=length( V10))
colnames(df.read_per_txpeak)=c("peak_chrom", "peak_start", "peak_end", "dataset", "annot_transcript_id", "peakid", "readcount_link")
write.table(df.read_per_txpeak, args[2], quote=F, sep="\t", col.names = T, row.names = F)