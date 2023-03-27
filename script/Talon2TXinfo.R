args=commandArgs( trailingOnly=T )
library(dplyr)
library(tidyr)
options(scipen=999)

df=read.table(args[1], stringsAsFactors=F, header=T)
mincount=as.numeric(args[2])

df.no_read_info=df[, !colnames(df) %in% c("read_name", "chrom", "read_start", "read_end", "genome_build", "read_length")]
df.no_read_info<-df.no_read_info %>% group_by( dataset, annot_transcript_id  ) %>% dplyr::mutate( readcount=length(annot_transcript_id));
df.no_read_info<-df.no_read_info %>% group_by( dataset, annot_transcript_id ) %>% dplyr::slice_sample(n=1)
df.no_read_info=df.no_read_info[,c("annot_transcript_id", "annot_gene_id", "gene_novelty" , "transcript_novelty", "ISM_subtype", "annot_gene_name", "annot_transcript_name", "n_exons", "gene_ID", "transcript_ID", "dataset", "readcount")]



vecTX.mincount=unique( df.no_read_info[df.no_read_info$readcount>=mincount, ]$annot_transcript_id )


#df.libsize=read.table(args[2], stringsAsFactors=F, header=T)
#dfm=merge( df.no_read_info, df.libsize[, c("dataset", "libsize")], by="dataset")
#dfm$TPM=round( (dfm$readcount*10^6)/dfm$libsize, 2)
#dfmg<-df.no_read_info %>% tidyr::spread( dataset, readcount, fill=0)
df.no_read_info$PassMinCount=ifelse(df.no_read_info$annot_transcript_id %in% vecTX.mincount, "O", "X")


write.table( df.no_read_info, args[3], col.names=T, row.names=F, sep="\t", quote=F)
 



