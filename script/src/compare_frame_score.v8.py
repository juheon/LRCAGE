import math
import sys
import re
#from Bio import SeqIO
import pysam
from Bio import Align
from Bio.Align import substitution_matrices


#For each read, compare its CDS against all CDS isoforms of Ref genes.
## as a result, one read may have multiple frame call. then choose as follows (normal > truncated > chimeric-normal > chimeric-truncated)
ENUM_NORMAL=1
ENUM_TRUNC=2
ENUM_CHIMNORM=3
ENUM_CHIMTRUN=4
ENUM_OUT=5
ENUM_NEW=6
ENUM_NULL=7

SCORE_CUTOFF=100	#approximately 15 aminoacids match
FRAC_MATCH=0.5		#Min fraction of matched bases out of entire bases in query sequence

thread=16

bByTALON=True

bTest=False
vecTestTX=["PB.10823_prot_1"]

icds_query=sys.argv[1]	#PB protein id, PB transcirpt id
icds_ref=sys.argv[3] 	#ref transcript id
imap=sys.argv[2] 	#PB transcript id, ref gene
oprefix=sys.argv[4] 

#Find the corresponding GENCODE genes using TALON
#prot to ref gene id
dictProt2RefGene={}
bHeader=True
nIndex_txname=nIndex_gencodeid=nIndex_protein=-1
for line in open(imap):
	arr=line.rstrip("\n").split("\t")
	if bHeader:
		bHeader=False;
		continue	
	
	protid=arr[2]
	geneid=arr[0]

	dictProt2RefGene[ protid ]=geneid
	
print("Load Prot to Ref gene ID")


def ExtractTXID(a_strname):
	arr=a_strname.split("|")
	return(arr[1]);	

def ExtractGeneID(a_strname):
	arr=a_strname.split("|")
	return(arr[2])


#Load indexes of Ref CDS
##For each gene, find the longest CDS
dictRef={} #nested dictionary
for record in pysam.FastxFile( icds_ref ):
	strTXID=record.name;
	strGeneID=record.comment;	
	if strGeneID not in dictRef:
		dictRef[ strGeneID ]=[]
		dictRef[ strGeneID ].append( (strTXID, record ));
	else:
		dictRef[ strGeneID ].append( (strTXID, record ));
print("Load Ref CDS index")

def ConvertENUM2string(a_enum):
	if a_enum==ENUM_NORMAL:
		return "NORMAL"
	elif a_enum==ENUM_TRUNC:
		return "TRUNCATED"
	elif a_enum==ENUM_CHIMNORM:
		return "CHIMERIC_NORMAL";
	elif a_enum==ENUM_CHIMTRUN:
		return "CHIMERIC_TRUNCATED";
	elif a_enum==ENUM_NEW:
		return "UNANNOTATED_CDS";
	elif a_enum==ENUM_OUT:
		return "OUTFRAME"
	else:
		return "NA"

##Use the default value of EMBOSS Water for protein alignment
##EMBOSS Water algorithm uses BLOSUM62 which 
##higher BLOSUM nubmer is for more closely related proteins
##BLOSUM62 is built from >=62% identitical proteins (measure the frequency of amino acid change)
##I will use BLOSUM90 to reduce the spurious alignment with multipel .... (substitution)


aligner = Align.PairwiseAligner()
aligner.mode="local"    #Smith-Waterman
matrix = substitution_matrices.load("BLOSUM90")
matrix = matrix.select(matrix.alphabet + "U")
aligner.substitution_matrix = matrix	#https://github.com/biopython/biopython/issues/3205

aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def ConvertAlignCoord2AACoord( a_arrIdx, a_idx_first, a_dictDel, a_range=False):	
	if len(a_arrIdx)==0:
		return [];

	if a_range==False:
		arrNew=[];
		for i in a_arrIdx:
			arrNew.append( i-a_idx_first-( a_dictDel[ i ] if i in a_dictDel else 0) );
		return( arrNew )
	else:
		arrNewRanges=[];
		for i in a_arrIdx:
			newRange=[ i[0]-a_idx_first-( a_dictDel[ i[0] ] if i[0] in a_dictDel else 0), i[1]-a_idx_first-( a_dictDel[ i[1] ] if i[1] in a_dictDel else 0) ]
			arrNewRanges.append( newRange );
		return( arrNewRanges );


def GetFrameCall( a_qrecord, a_refrecord ):
	eType5end=eType3end=ENUM_NULL;
	align_final=None;

	strQseq=a_qrecord.sequence.rstrip("*");
	nLen_query=len( strQseq );

	if a_refrecord is None:
		##If ref gene is non-coding RNA
		return (ENUM_NEW, ENUM_NEW, None, nLen_query, -1, -1, -1, -1, -1, -1, -1, [] );

	try:
		strRseq=a_refrecord.sequence;
		nLen_ref=len( strRseq )
		alignments=aligner.align( strQseq, strRseq )

		if len(alignments)<1:
			#No alignment 
			eType5end=eType3end=ENUM_OUT
			return (eType5end, eType3end, None, nLen_query, nLen_ref, -1, -1, -1, -1, -1, -1, []);

		align_final=alignments[0];
		strQueryPerBase=str(align_final).split("\n")[0]
		strAlignPerBase=str(align_final).split("\n")[1]
		strSubjectPerBase=str(align_final).split("\n")[2]

		##Find indexes of gap, mismatch relative to the first alinged AAs (later this will be transformed into Query AA coordinate)
		##	MVER		query
		##	||||		align	
		##MHARGEMVERAGTHATH	subject
		##
		arrIdx_align_inalign=[x.start() for x in re.finditer("-|\\.|\\|", strAlignPerBase)];
		nIdx_align_first=min( arrIdx_align_inalign );
		nIdx_align_last=max( arrIdx_align_inalign );

		#arrIdx_space_inquery=[x.start() for x in re.finditer(" ", strQueryPerBase)]
		#arrIdx_space_inalign=[x.start() for x in re.finditer(" ", strAlignPerBase)]
		
		##mismatch
		arrIdx_mismatch=[x.start() for x in re.finditer("\\.", strAlignPerBase)]

		##insertion deletion
		arrIdx_gap_inquery=[x.start() for x in re.finditer("-", strQueryPerBase)]
		arrIdx_gap_inalign=[x.start() for x in re.finditer("-", strAlignPerBase)]
		arrIdx_del=list( set( arrIdx_gap_inquery ) & set( arrIdx_gap_inalign ) ); arrIdx_del.sort();
		arrIdx_ins=list( set( arrIdx_gap_inalign ) - set( arrIdx_gap_inquery ) ); arrIdx_ins.sort();

		##chimeric
		arrIdx_letter_inquery=[x.start() for x in re.finditer("[A-z]", strQueryPerBase )]
		
		arrIdx_chim=list( set(arrIdx_letter_inquery) - set(arrIdx_align_inalign) ); arrIdx_chim.sort()
		arrIdx_chim5=[ x for x in arrIdx_chim if x<nIdx_align_first ] if len(arrIdx_chim)>0 else [];
		arrIdx_chim3=[ x for x in arrIdx_chim if x > nIdx_align_last ] if len(arrIdx_chim)>0 else [];

		##truncated
		arrIdx_letter_insubject=[x.start() for x in re.finditer("[A-z]", strSubjectPerBase)];
		arrIdx_trun=list( set(arrIdx_letter_insubject) - set(arrIdx_align_inalign) )

		arrIdx_trun5=[];
		arrIdx_trun3=[];
		if len( arrIdx_trun )>0:
			arr=[ x for x in arrIdx_trun if x<nIdx_align_first ]		
			arrIdx_trun5=[ max( arr) ] if len(arr)>0 else [];

			arr2=[ x for x in arrIdx_trun if x>nIdx_align_last ]
			arrIdx_trun3=[ min(arr2) ] if len(arr2)>0 else [];		

		##For insertion, one indel is defined as range
		arrIdx_ins_group=[];
		if len(arrIdx_ins)>0:
			arrIdx_notconsecutive=[];
			arrIdx_notconsecutive=[ i for i in range(0, len( arrIdx_ins )-1) if arrIdx_ins[i+1]-arrIdx_ins[i]>1] + [len(arrIdx_ins)-1];
			for i in range( 0, len(arrIdx_notconsecutive) ):
				if i==0:
					arrIdx_ins_group.append( [ arrIdx_ins[ 0 ], arrIdx_ins[ arrIdx_notconsecutive[i] ] ]);
				else:
					arrIdx_ins_group.append( [ arrIdx_ins[ arrIdx_notconsecutive[i-1]+1 ], arrIdx_ins[ arrIdx_notconsecutive[i] ] ]);
		
		arrIdx_del_group=[];
		if len(arrIdx_del)>0:
			arrIdx_notconsecutive=[];
			arrIdx_notconsecutive=[ i for i in range(0, len( arrIdx_del )-1) if arrIdx_del[i+1]-arrIdx_del[i]>1] + [len(arrIdx_del)-1];
			for i in range( 0, len(arrIdx_notconsecutive) ):
				if i==0:
					arrIdx_del_group.append( [ arrIdx_del[ 0 ], arrIdx_del[ arrIdx_notconsecutive[i] ] ]);
				else:
					arrIdx_del_group.append( [ arrIdx_del[ arrIdx_notconsecutive[i-1]+1 ], arrIdx_del[ arrIdx_notconsecutive[i] ] ]);
		

		##Dictionary, aggregated length of deletions before the align coordinate
		dictAggDel={};
		if len(arrIdx_del)>0:
			for i in range( min(arrIdx_letter_inquery), max(arrIdx_letter_inquery)+2 ):#Extend one extra bp for 3' truncated  
				arrAggDel_new=list( set( range(0, i) ).intersection( set(arrIdx_del) ) )
				arrAggDel_new.sort();
				dictAggDel[ i ]=len( arrAggDel_new )
		
			


		#Convert align coordinate into query coordinate
		##If there is deletion, align coordinate to protein coordinate should be accounted for deletion size
                ##      MVE---T            query
                ##      |||---T            align
                ##MHARGEMVERAGTHATH        subject
		##0123456789X		   align coord
		##	      *		align coord: 12, query coord: 3
		##			12-6-3=3

		#All bed of VAR region is defined as the exact location of VAR
		##However, for TRUNCATED, the location is VAR is coordiate outside the VAR region
		##This lead to not detection by bedtools intersect command in the downstream
		##For this reason, I will use the first and the last coordinate of protein as the TRUCNATED position

		nIdx_query_start=min(arrIdx_letter_inquery)
		arrQIdx_mismatch=ConvertAlignCoord2AACoord( arrIdx_mismatch, nIdx_query_start, dictAggDel );
		arrQIdx_ins_group=ConvertAlignCoord2AACoord( arrIdx_ins_group, nIdx_query_start, dictAggDel, True );
		arrQIdx_del_group=ConvertAlignCoord2AACoord( arrIdx_del_group, nIdx_query_start, dictAggDel, True );
		arrQIdx_chim5=ConvertAlignCoord2AACoord( arrIdx_chim5, nIdx_query_start, dictAggDel );
		arrQIdx_chim3=ConvertAlignCoord2AACoord( arrIdx_chim3, nIdx_query_start, dictAggDel );
		arrQIdx_trun5=ConvertAlignCoord2AACoord( arrIdx_trun5, nIdx_query_start, dictAggDel );
		arrQIdx_trun3=ConvertAlignCoord2AACoord( arrIdx_trun3, nIdx_query_start, dictAggDel );
	
		arrQIdx_combined=[];
		for i in arrQIdx_mismatch:
			arrQIdx_combined.append( [i,i+1, "MISMATCH", 1] );
		for i in arrQIdx_ins_group:
			arrQIdx_combined.append( [i[0],i[1]+1, "INSERTION", i[1]+1-i[0] ] );
		for i in arrQIdx_del_group:
			arrQIdx_combined.append( [i[0]-1,i[1]+1, "DELETION", 2 ] );
		if len(arrQIdx_chim5)>0:
			arrQIdx_combined.append( [arrQIdx_chim5[0], arrQIdx_chim5[-1]+1, "CHIMERIC_5", arrQIdx_chim5[-1]+1-arrQIdx_chim5[0] ] );
		if len(arrQIdx_chim3)>0:
			arrQIdx_combined.append( [arrQIdx_chim3[0], arrQIdx_chim3[-1]+1, "CHIMERIC_3", arrQIdx_chim3[-1]+1-arrQIdx_chim3[0] ] );
		if len(arrQIdx_trun5)>0:
			arrQIdx_combined.append( [arrQIdx_trun5[0]+1, arrQIdx_trun5[0]+2, "TRUNCATED_5", 1] );
		if len(arrQIdx_trun3)>0:	
			arrQIdx_combined.append( [arrQIdx_trun3[0]-1, arrQIdx_trun3[0], "TRUNCATED_3", 1] );

		#old version	
		##I don't replace these code with newer code using coordinate information. but I will just make sure the information from old version is the same with the new version
		nLen_align=len(strAlignPerBase.strip(" "))
		nLen_align_match=strAlignPerBase.count("|")
		nLen_align_gap=strAlignPerBase.count("-")
		nLen_align_mismatch=strAlignPerBase.count(".");
		nLen_query_novel5=0
		nLen_query_novel3=0
		fMatch=nLen_align_match/nLen_align

		if bTest:
			print(align_final)
			print(align_final.score)


		arrStatus_query=align_final.aligned[0]
		arrStatus_ref=align_final.aligned[1]

		#5end, 3end relative to genome (strand insensitive)
		##5end: Does leftmost base of query is mapped to reference
		bInclude5end_query=(arrStatus_query[0][0]==0)           #>0: 5' end of query is novel
		bInclude3end_query=(arrStatus_query[-1][1]==nLen_query)#<0: 3' end of query is novel
	
		bInclude5end_ref=(arrStatus_ref[0][0]==0)               #>0: 5' end of ref is not covered by query
		bInclude3end_ref=(arrStatus_ref[-1][1]==nLen_ref)       #<0: 3' end of ref is not covered by query
	
		#nInclude5end_query     nInclude5end_ref        description
		#T                      T                       perfect match
		#T                      F                       left truncated
		#F                      T                       chimeirc normal
		#F                      F                       chimeric truncated

		if not bInclude5end_query:
			nLen_query_novel5=arrStatus_query[0][0]
		if not bInclude3end_query:
			nLen_query_novel3=nLen_query-arrStatus_query[-1][1]

		#test code
		if bTest:
			if nLen_align_mismatch!=len(arrQIdx_mismatch):
				print("mismatch error"+"\t"+a_qrecord.name);
			if (nLen_query_novel5>0) != (len(arrQIdx_chim5)>0):
				print("chim5 error"+"\t"+a_qrecord.name);
			if (nLen_query_novel3>0) != (len(arrQIdx_chim3)>0):
				print("chim3 error"+"\t"+a_qrecord.name);
			if ((len(arrQIdx_ins_group)>0) | (len(arrQIdx_del_group)>0) ) != (nLen_align_gap>0):
				print("gap error"+"\t"+a_qrecord.name+"\t"+str(len(arrQIdx_ins_group))+"\t"+ str(len(arrQIdx_del_group))+"\t" +str(nLen_align_gap));		


		if (fMatch < FRAC_MATCH) or (align_final.score < SCORE_CUTOFF): 
			#too small match
			eType5end=eType3end=ENUM_OUT
			return (eType5end, eType3end, None, nLen_query, nLen_ref, nLen_align_match, nLen_align_mismatch, nLen_align_gap, nLen_query_novel5, nLen_query_novel3, fMatch, []);

		if bInclude5end_query and bInclude5end_ref:
			eType5end=ENUM_NORMAL;	
		elif bInclude5end_query and not bInclude5end_ref:
			eType5end=ENUM_TRUNC
		elif not bInclude5end_query and bInclude5end_ref:
			eType5end=ENUM_CHIMNORM
		else:
			eType5end=ENUM_CHIMTRUN

		if bInclude3end_query and bInclude3end_ref:
			eType3end=ENUM_NORMAL;
		elif bInclude3end_query and not bInclude3end_ref:
			eType3end=ENUM_TRUNC
		elif not bInclude3end_query and bInclude3end_ref:
			eType3end=ENUM_CHIMNORM
		else:
			eType3end=ENUM_CHIMTRUN

		#5end		3end		
		#M		M		
		#M		T		
		#M		CN		
		#M		CT		
		#T		M		
		#T		T		
		#T		CN		
		#T		CT		
		#CN		M		
		#CN		T 		
		#CN		CN		
		#CN		CT		
		#CT		M		
		#CT		T
		#CT		CN		
		#CT		CT		
		
		return (eType5end, eType3end, align_final, nLen_query, nLen_ref, nLen_align_match, nLen_align_mismatch, nLen_align_gap, nLen_query_novel5, nLen_query_novel3, fMatch, arrQIdx_combined);
	except OverflowError:
		eType5end=eType3end=ENUM_NULL
		return (eType5end, eType3end, None, nLen_query, nLen_ref, -1, -1, -1, -1, -1, -1, []);	
#	except ValueError:
#		print("ERR "+a_qrecord.id+"\t"+a_qrecord.seq)
#		return (eType5end, align_final);

def ConvertScores2String(a_qlen, a_rlen, a_qmatch, a_qmismatch, a_qgap, a_novel5, a_novel3):
	#query_len(match, mismatch, gap)

	return str(a_qlen)+" (M:"+str(a_qmatch)+",MM:"+str(a_qmismatch)+",G:"+str(a_qgap)+",N5:"+ str(a_novel5)+",N3:"+str(a_novel3)+");"+str(a_rlen);


#Load all read
listRead=[];
for record in pysam.FastxFile( icds_query  ):
	strProtName=record.name;	
	strTXName=record.comment

	if bTest and strProtName not in vecTestTX:
		continue;

	listRead.append( (strProtName, strTXName, record) )


#Split read into multiple group for parallel processing
nTotal=len(listRead)
nGroup=nSizePerGroup=0;
if nTotal % thread == 0:
	nSizePerGroup=nTotal/thread;
	nGroup=int(nTotal/nSizePerGroup);
else:
	nSizePerGroup=int( math.ceil( nTotal/thread ));
	nGroup=int( math.ceil(nTotal/nSizePerGroup));


listReadgroup=[];
for i in range(0, nGroup ):
	nStartIdx=int( nSizePerGroup*i )
	nEndIdx=int( nSizePerGroup*i+nSizePerGroup if i<nGroup else nTotal-1 )
	listReadgroup.append( listRead[nStartIdx:nEndIdx])
print("Prepare multiprocessing")

def GetArrTXAgainst( a_txquery, a_dictProt2RefGene, a_dictRef ):
	#Find the list of annotated protein sequence for comparison	

	#find ref gene id corresponding to protein id		
	if a_txquery not in a_dictProt2RefGene:
		print("invalid query transcript " + a_txquery)
		return [];

	strRefGene=a_dictProt2RefGene[ a_txquery ]
	#find all ref transcripts associated with ref gene
	if strRefGene not in a_dictRef:
		print("invalid gene " + strRefGene)
		return [];

	arrRefRecord=a_dictRef[ strRefGene ]
	return arrRefRecord

def GetRepFrameCall( a_queryrecord, a_refrecord ):
	tpFrame_final=None
	strTX_final=""
	if a_refrecord is None:
		#no CDS that compared against; UNANNOCATED_CDS
		tpFrame_final=GetFrameCall( a_queryrecord, None )
		strTX_final="NA"
	else:
		#find the best matching canddiate
		eFrame5_final=eFrame3_final=ENUM_NULL;
		tpFrame_final=None; nScore_final=-1;

		for tuplRecord in a_refrecord:
			reftxid=tuplRecord[0]; refrecord=tuplRecord[1];
			tpFrame=GetFrameCall( a_queryrecord, refrecord );
			eFrame5=tpFrame[0]; eFrame3=tpFrame[1];
			nScore=0 if tpFrame[2] is None else tpFrame[2].score;
			if nScore > nScore_final:
				eFrame5_final=eFrame5; eFrame3_final=eFrame3;
				tpFrame_final=tpFrame;  nScore_final=nScore;
				strTX_final=reftxid;
	if tpFrame_final is None:
		print("impossible "+strProtName);
	return (tpFrame_final, strTX_final);		


def FindFinalFrameCall( a_listRead ):
	arrFinalFrame=[];	#tpFrame, 

	for i in a_listRead:
		strProtName=i[0]
		strTXName=i[1]
		record=i[2]	


		#strTXname=record.comment.split(";")[0]	
		#Find corresponding RefGene of protein
		strRefGene=dictProt2RefGene[ strProtName ] if strProtName in dictProt2RefGene else None;
		if strRefGene is None:
			print("ERR")
			continue;	
	
		#Find all transcripts of RefGene
		arrRecordRefTx=dictRef[ strRefGene ] if strRefGene in dictRef else None;

		tpFrame_final=GetRepFrameCall(record, arrRecordRefTx)
		arrFinalFrame.append( (strProtName, strTXName, len(record.sequence), tpFrame_final[0], tpFrame_final[1] ));
	return( arrFinalFrame )	
	





#multithreading
from multiprocessing import Pool
p = Pool( thread )
result=p.map( FindFinalFrameCall, listReadgroup );
p.close();
p.join();

print("ORF comparison")



#print to files
#ofile=open(oprefix+".framecall.txt", "w+")
ofile_align=open(oprefix+".framecall_align.txt", "w+")
#ofile_coord=open(oprefix+".framecall_diff.bed", "w+");  #coordinates of chimeric, truncated, outframe, gap, mismatch

results=[ j for i in result for j in i ]; ##flatten arrays
#ofile.write( "\t".join( ["readname", "class", "class_3end", "ref_geneID", "ref_txID", "pct_match", "len_read", "len_reftx", "match", "mismatch", "gap", "novel5", "novel3"]  )+"\n" );
for i in results:
	strReadName=i[0];	
	strReadTX=i[1]
	nSeqLen=i[2];
	tpFrame=i[3]
	strTX=i[4]

	varRegions=[];

	eFrame5=eFrame3=ENUM_NULL
	fMatch=0; nLen_query=str(nSeqLen); nLen_ref=0;
	nLen_align_match=nLen_align_mismatch=nLen_align_gap=nLen_query_novel5=nLen_query_novel3=0;
	alignPlot="";
	
	if tpFrame is None:
		print("ERROR tpFrame is none "+strReadName );
		continue;


	eFrame5=tpFrame[ 0 ]; eFrame3=tpFrame[ 1 ] 
	
	fMatch=tpFrame[10]	
	nLen_query=tpFrame[3]; 	nLen_ref=tpFrame[4]
	nLen_align_match=tpFrame[5]; 	nLen_align_mismatch=tpFrame[6]
	nLen_align_gap=tpFrame[7]
	nLen_query_novel5=tpFrame[8]; 	nLen_query_novel3=tpFrame[9]

	alignPlot=tpFrame[2]

	varRegions=tpFrame[ 11 ]


	arrAlignInfo=[ strReadName, strReadTX, strTX, "LenQuery: "+str(nLen_query), "LenRef: "+str(nLen_ref) ]

#	print(arrAlignInfo)
	#ofile.write( "\t".join( arrAlignInfo )+"\n" );


	#align infomration
	strAlignDetail=str(alignPlot);
	if strAlignDetail=="None":
		strAlignDetail="empty\nempty\nempty\n"

	ofile_align.write(",".join( arrAlignInfo )+"\n"+strAlignDetail+"\n");		
	
	#variable region	
##	if len( varRegions )<1:
##		#either OUTFRAME or NEW CDS
##		strVar=ConvertENUM2string( eFrame5 );
##		ofile_coord.write("\t".join([ strReadName, str(0), str(nSeqLen), strVar, str(nSeqLen), "+" ])+"\n" );
##
##	for j in varRegions:	
##		ofile_coord.write("\t".join([ strReadName, str(j[0]), str(j[1]), j[2], str(j[3]), "+" ])+"\n" );
print("Print to files")

#ofile.close();
ofile_align.close();
#ofile_coord.close();




print("done");





