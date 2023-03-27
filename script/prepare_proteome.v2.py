import pysam
import sys
import copy

ipep=sys.argv[1]
italon=sys.argv[2]
oprefix=sys.argv[3]

#if the multiple isoform of of the same produces the same protein sequence, they are combined
#output
# combined pep
# transcript id to protein id mapping table by TALON

#
bFirst=True;
dictTX2Gene={}
nIdxGeneId=nIdxTXId=-1
for line in open(italon):
	arr=line.rstrip("\n").split("\t")
	if bFirst==True:
		bFirst=False;
		nIdxGeneId=arr.index("annot_gene_id")
		nIdxTXId=arr.index("annot_transcript_id")
		continue;
	
	txid=arr[nIdxTXId]
	if txid not in dictTX2Gene:
		dictTX2Gene[ txid ]=arr[ nIdxGeneId ]

		

#protein of "ENSG00000044115.20_prot_1" is missing in the output why?
with pysam.FastxFile( ipep ) as fin, open( oprefix+".proteome.fasta" , "w") as fout, open(oprefix+".proteome_table.txt", "w+") as fout2:
	dictGeneTx={}

	#Load all protein sequence 
	for entry in fin:
		strTxname=entry.name.split(":")[0]
		print(entry.name)

		strGenename=dictTX2Gene[ strTxname ];		
		if strGenename not in dictGeneTx:
			dictGeneTx[ strGenename ]=[]

		dictGeneTx[ strGenename ].append(entry)

	fout2.write("\t".join(["#genename", "txname", "protName"])+"\n")

	dictTx2ProtTable={}
	for genename in dictGeneTx:		
		#Sort transcript by protein size
		dictGeneTx[ genename ].sort( key=lambda x:len(x.sequence), reverse=True);	

		#Across all transcirpt in one gene, see if it has the same sequence as other transcript
		## If it is assign protein ID.
		dictTx2ProtTable[ genename ]={}
		dictReportedSeq={}	#seq,  (protein name, [transcript id])
		for entry2 in dictGeneTx[ genename ]:
			strTxname=entry2.name.split(":")[0];
			strProtName="NA"
			if entry2.sequence in dictReportedSeq:
				strProtName=dictReportedSeq[ entry2.sequence ][ 0 ]
				dictReportedSeq[ entry2.sequence ][1].append( strTxname )
			else:
				strProtName=genename+"_prot_"+str( len( dictReportedSeq )+1)
				dictReportedSeq[ entry2.sequence ]=( strProtName, [strTxname] );
			dictTx2ProtTable[ genename ][ strTxname ]=strProtName		
			fout2.write("\t".join([genename, strTxname, strProtName])+"\n")	
		
		##For each sequence, we concat all transcripts and output as fasta header
		dictOutProtein={}
		for entry2 in dictGeneTx[ genename ]:
			strTxname=entry2.name.split(":")[0];
			strProt=dictReportedSeq[ entry2.sequence ][0]
				
			if strProt not in dictOutProtein: #given protein is not outputted to fasta file yet
				entry2.name=strProt;
				entry2.comment=";".join(dictReportedSeq[ entry2.sequence ][1])
				fout.write(str(entry2)+"\n");

				dictOutProtein[ strProt ]=strProt
			else:
				continue;






			



