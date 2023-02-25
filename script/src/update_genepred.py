
import sys 
#import pysam


igenepred=sys.argv[1]
txinfo=sys.argv[2]
oprefix=sys.argv[3]
bPeakExtend=sys.argv[4]	#1: True, 0: False

#TALON genepred extend/shorten 1st exon if the assembled transcript match well with GENCODE annotation
#If peak, it override the gene annotation
#If the same transcript in two different dataset has two different peaks, the upper most one is chosen as the representative

dictPeakPerConfTx={}
dictPeakPerConfTx_New={}
nIndexTxId=-1; nIndexPeakId=-1; 
nIndexDataset="";nIndexPeakChrom=""; 
nIndexPeakStart=-1;nIndexPeakEnd=-1; nIndexTXnovelty=-1
bFirst=True;

print("Load Transcript Info")
nIndex=0
nIndex_IsConfident=-1;
for line in open( txinfo ):
	arr=line.rstrip("\n").split("\t")

	if bFirst==True:
		nIndexPeakId=arr.index( "peakid" )
		nIndexTxId=arr.index("annot_transcript_id")
		nIndexPeakChrom=arr.index("peak_chrom")
		nIndexPeakStart=arr.index("peak_start")
		nIndexPeakEnd=arr.index("peak_end")
		nIndexTXnovelty=arr.index("transcript_novelty");
		nIndex_IsConfident=arr.index("IsConfident")
		bFirst=False
		continue;	

	strIsConfident=arr[ nIndex_IsConfident ]
	if strIsConfident=="X":
		continue;

	txid=arr[ nIndexTxId ]
	peakid=arr[ nIndexPeakId ]
	arrpeakcoord=[ peakid, arr[nIndexPeakChrom], int(arr[nIndexPeakStart]), int(arr[nIndexPeakEnd]) ]


	if txid not in dictPeakPerConfTx:
		dictPeakPerConfTx[txid]=arrpeakcoord	

	if arr[nIndexTXnovelty]=="Known":
		continue;
	if txid not in dictPeakPerConfTx_New:
		dictPeakPerConfTx_New[txid]="";


#Store the number of exons per transcript for reverse stranded transcript
def GetFirstExon(a_strstart, a_strend, a_strand):
	arrstart=a_strstart.rstrip(",").split(",");
	arrend=a_strend.rstrip(",").split(",");

	if a_strand=="+":
		return (int( arrstart[0]), int(arrend[0]));
	else:
		return (int( arrstart[-1]), int(arrend[-1]))


def IsPeakInExon(a_peakstart, a_peakend, a_exonstart, a_exonend):
	setPeak=set( range( a_peakstart+1, a_peakend+1 ) );
	setExon=set( range( a_exonstart+1, a_exonend+1 ) );

	nPeakSize=a_peakend-a_peakstart
	return True if len( setPeak.intersection( setExon ) )==nPeakSize else False;

print("Update genepred")
ifile=open( igenepred )
ofile=open(oprefix+".conf.genepred", "w+")
ofile2=open(oprefix+".conf_nc.genepred", "w+")
dictLastExon={}
for entry in ifile:
	arr=entry.rstrip("\n").split("\t")
	txid=arr[0]

	if txid not in dictPeakPerConfTx:	
		#print("missing peak")
		continue;

	strStrand=arr[2]	

	if bPeakExtend=="True":	
		#Exon boundary update
		arrpeak=dictPeakPerConfTx[ txid ];
		tpFirstExon=GetFirstExon( arr[8], arr[9], strStrand );

		nStartPeak=arrpeak[2]
		nEndPeak=arrpeak[3]
	
		bIsPeakInExon=IsPeakInExon( nStartPeak, nEndPeak, tpFirstExon[0], tpFirstExon[1])
		if bIsPeakInExon:
			#if entire peak is located within exon, the exon is shorten
			if strStrand=="+":
				nStart=max( int( arr[3]), nStartPeak ); arr[3]=str(nStart);
				brr=arr[8].split(",");
				brr[0]=str(nStart);
				arr[8]=",".join(brr);
			else:
				nEnd=min( int(arr[4]), nEndPeak ); arr[4]=str(nEnd);
				brr=arr[9].split(",");
				brr[-2]=str(nEnd);
				arr[9]=",".join(brr);
		else:
			#if some of peak is located out of exon, the exon is extended
			if strStrand=="+":
				nStart=min( int(arr[3]), nStartPeak ); arr[3]=str(nStart);
				brr=arr[8].split(",");
				brr[0]=str(nStart);
				arr[8]=",".join(brr);
			else:
				nEnd=max( int(arr[4]), nEndPeak ); arr[4]=str(nEnd);
				brr=arr[9].split(",");
				brr[-2]=str(nEnd);
				arr[9]=",".join(brr);
	
		arr[5]=arr[6]=arr[4];

	ofile.write("\t".join(arr)+"\n");

	if txid in dictPeakPerConfTx_New:
		ofile2.write("\t".join(arr)+"\n");

ofile.close();ofile2.close();
ifile.close()







