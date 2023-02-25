args=commandArgs(trailingOnly=TRUE)
library(dplyr)

inoveltx=args[1]
itx2prot=args[2]   #which tx become which protein
iprot2type=args[3] #is protein chimeric or normal against GENCODE proteins
output=args[4]
df.txinfo=read.table(inoveltx, stringsAsFactors = F, header=T)


#For each gene, what are transcripts and their proteinsll -
df.tx2prot=read.table(itx2prot, stringsAsFactors = F, header=F);
df.tx2prot$txid=apply(df.tx2prot, 1, function(a_x){
  strTX=stringr::str_split(a_x[2], "\\(")[[1]][1]
  return( strTX )
})
colnames(df.tx2prot)[1:3]=c("geneid", "txid_long", "protid")

#At the level of protein, their comparison to GENCODE proteins
df.stat=read.table(iprot2type, stringsAsFactors = F, header=T)

dfm=merge(df.tx2prot[, c("protid", "txid", "geneid")], df.stat, 
          by.x="protid", by.y="readname", all=T)


dfmt<-merge( df.txinfo, dfm, by.x="annot_transcript_id", by.y="txid", all.x=T)

write.table(dfmt, output, sep="\t", quote=F, col.names = T, row.names = F)


