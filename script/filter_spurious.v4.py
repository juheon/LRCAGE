import sys
import pysam
import math
bFirst=True
dictAllTX={}
dictPair={}
ifalign= open(sys.argv[1]);
query_gtf=sys.argv[2]
ref_gtf=sys.argv[3]
ffile=open(sys.argv[4], "w+");

dictQueryGTF={}
dictRefGTF={}

#1. Load align data (protein name, transcript name)
#2. Load CDS coordinate of query gtf and ref gtf
#3. For each protein in 1, find the transcript in query gtf
#4. 

#query gtf, ref gtf to find CDS coordinate


bTest=False
arrTest=["PB.6509.52", "PB.1368.23"]
arrTest=["PB.10823.17"]
arrTest=["PB.7637.12"]
arrTest=["PB.5212.16"]
arrTest=["TALONT000261326"]
arrTest=["TALONT000300662"]
with pysam.TabixFile( query_gtf ) as fqgtf, pysam.TabixFile( ref_gtf ) as frgtf:
	print("Load Query GTF")
	for entry in fqgtf.fetch( parser=pysam.asGTF() ):
		if entry.feature!="CDS":
			continue;

		if bTest and entry.transcript_id not in arrTest:
			continue;
		
		if entry.transcript_id not in dictQueryGTF:
			dictQueryGTF[ entry.transcript_id ]=[entry.strand, [i+1 for i in range(entry.start, entry.end )]]
		else:
			dictQueryGTF[ entry.transcript_id ][1]=dictQueryGTF[ entry.transcript_id ][1]+[i+1 for i in range(entry.start, entry.end )]

	print("Load Ref GTF")
	for entry in frgtf.fetch( parser=pysam.asGTF() ):
		if entry.feature!="CDS":
			continue;

#		if bTest and entry.transcript_id not in ["ENST00000392590.3", "ENST00000271688.10", "ENST00000358389.7", "ENST00000543610.5", "ENST00000288757.7", "ENST00000534624.5"]:
#			continue
	
		if entry.transcript_id not in dictRefGTF:
			dictRefGTF[ entry.transcript_id ]=[entry.strand, [i+1 for i in range(entry.start, entry.end )]]
		else:
			dictRefGTF[ entry.transcript_id ][1]=dictRefGTF[ entry.transcript_id ][1]+[i+1 for i in range(entry.start, entry.end )]

	for txid in dictQueryGTF:
		dictQueryGTF[ txid ][1].sort();
	for txid in dictRefGTF:
		dictRefGTF[ txid ][1].sort();



def GetArrGcoordByAApos(a_txid, a_dic, a_npos):

	#For given AA pos, return 3 gcoords	
	if a_txid not in a_dic:
		print("[GetArrGcoordByAApos] Not found: "+a_txid)
		return None;

	strStrand=a_dic[ a_txid ][0]
	arrGcoord=a_dic[ a_txid ][1]
	nCDSexon=len( arrGcoord )

	nLenProt=int( nCDSexon/3 )
	nGcoor_start=-1;nGcoor_end=-1;
	if strStrand=="+":
		nGcoor_start=3*(a_npos-1)+1
		nGcoor_end=3*(a_npos-1)+3
	else:
		nGcoor_start=3*(nLenProt-a_npos)+1
		nGcoor_end=3*(nLenProt-a_npos)+3

	#if a_npos==432:
	#	print( "GetArrGcoordByAApos" )
	#	print( (a_txid, nLenProt, nCDSexon) )
	#	print( (nGcoor_start, nGcoor_end) )	
	#	print(  arrGcoord[ (nGcoor_start-1):nGcoor_end ] )
	
	return( arrGcoord[ (nGcoor_start-1):nGcoor_end ] )		

	




nMinMatch=1	#if -.-, -|-, replaced into ---
fCutoff=0.1
##def GetFiltered(a_stralign):
def Convert2chunk_byblank(a_stralign):
	if a_stralign=="empty":
		return [];

	if len(a_stralign)==0:
		return [];

	arrExon=[]      #list of arrBase
	arrPrevBase=[]
	nPrevMode=nCurMode=-1   #0: empty, 1: match, 2: -
	for i in a_stralign:
		if i==" ":
			nCurMode=0;
		else:
			nCurMode=1;
		if nPrevMode!=nCurMode and nPrevMode>=0:
			arrExon.append( arrPrevBase )
			arrPrevBase=[];
		arrPrevBase.append( i )
		nPrevMode=nCurMode
	arrExon.append( arrPrevBase )
	return arrExon;


def Convert2chunkarr( a_stralign, a_bignoreX=False ):
	if a_stralign=="empty":
		return [];

	arrExon=[]	#list of arrBase
	arrPrevBase=[]
	nPrevMode=nCurMode=-1	#0: empty, 1: match, 2: -
	for i in a_stralign:
		if i=="|" or i==".":
			nCurMode=1;
		elif i==" ":
			nCurMode=0;
		elif i=="-":
			nCurMode=2;
		elif i=="X":	#used for spurious exon removal
			if a_bignoreX:
				nCurMode=2;
			else:
				nCurMode=3;

		else:
			nCurMode=100
			print("ERR")

		if nPrevMode!=nCurMode and nPrevMode>=0:
			arrExon.append( arrPrevBase )
			arrPrevBase=[];

		arrPrevBase.append( i )
		nPrevMode=nCurMode	
	arrExon.append( arrPrevBase )
	return arrExon;







NUM_LARGE=10^10
def GetTrimmed( a_str, a_arrTrim, a_tuplAnchor):	
	#count number of "-" in the trimemd region at 5' end and 3' end
	arr_dash=[i for i in range(0, len(a_str)) if a_str[i]=="-"  ]
	arr_trim=[i for i in range(0, len(a_arrTrim)) if a_arrTrim[i]==True ]

	arr_dash_intrim=list(set( arr_dash).intersection(set(arr_trim)) )
	arr_dash_intrim.sort();
	nMinAnchor=a_tuplAnchor[0]
	nMaxAnchor=a_tuplAnchor[1]

	arr_5=[i for i in arr_dash_intrim if i<a_tuplAnchor[0]]
	arr_3=[i for i in arr_dash_intrim if i>a_tuplAnchor[1]]
	return( len(arr_5), len(arr_3) )	

def GetUpdateEnd( a_arrTrim, a_string, a_tupltrimlen ):
	nLen=len( a_string );
	nLen_trim=len( a_arrTrim )

	arrFilt=[]	
	for i in range(0, nLen):
		if i<nLen_trim:
			if a_arrTrim[i]==True and a_string[i]=="-":
				continue;
			arrFilt.append( a_string[i] );
		else:
			arrFilt.append( a_string[i] );
	arrFilt=[" "]*a_tupltrimlen[0]+arrFilt+[" "]*a_tupltrimlen[1]	
	return("".join(arrFilt))	

def list2bed(nums):
	#coldfix https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list/20725061#20725061
	nums = sorted(set(nums))
	gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
	edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
	return list(zip(edges, edges))


def slicebyindex( a_string, a_arrindex):
	nPrev=0;
	arrString=[]
	for i in a_arrindex:
		arrString.append( a_string[nPrev:(i+1)] );	
		nPrev=i+1
	if nPrev!=len(a_string):
		arrString.append( a_string[nPrev:len(a_string)] )
	return arrString
		

def insertbyindex( a_string, a_arrindex, a_arrinsert ):
	#add string in a_arrinsert at the index in a_arrindex 
	#one index after a_arrindex
	if len(a_arrindex)==0:
		return a_string;

	arrString=slicebyindex( a_string, a_arrindex )

	strResult=arrString[ 0 ];
	for i in range(0, len(a_arrinsert)):
		strResult=strResult+a_arrinsert[i]+arrString[i+1];

	return strResult

				

def RemoveSpuriousExon( a_strQuery, a_alignperbase, a_strSubject ):
	#Replace spurios exon into novel exon or skipped exon
#	arrIdx_markedspur=[ i for i in range(0, len( a_alignperbase)) if a_alignperbase[i]=="X" ];
#	if len(arrIdx_markedspur)==0:
#		return (a_strQuery, a_alignperbase, a_strSubject)
#	arrRange_markedspur=list2bed( arrIdx_markedspur );
#
#	arrIdx_dash=[ i for i in range(0, len( a_alignperbase)) if a_alignperbase[i]=="-" ];
#	arrRange_dash=list2bed( arrIdx_dash );


	#grouping spurious exons
	## -------XXXXX------XXXXX-----XXXXX----
	## if multipel spurious exons are located in intron, they are grouped into one
	arrIdx_query=[ i for i in range(0, len(a_strQuery)) if a_strQuery[i].isalpha() ]
	arrIdx_subject=[ i for i in range(0, len(a_strSubject)) if a_strSubject[i].isalpha() ]
	arrIdx_sus=[ i for i in range(0, len( a_alignperbase)) if a_alignperbase[i]=="X" ]
	arrIdx_chunk=[ i for i in range(0, len( a_alignperbase)) if a_alignperbase[i]=="X" or a_alignperbase[i]=="-" ] #Convert2chunkarr( a_alignperbase, True )
	arrRange_chunk=list2bed( arrIdx_chunk )
	if len(arrRange_chunk)==0:
		return (a_strQuery, a_alignperbase, a_strSubject)

	
	arrStrQ=list(a_strQuery); arrStrS=list(a_strSubject); arrStrA=list(a_alignperbase)
	arrRange_QSpur=[]; arrRange_SSpur=[]; arrRange_ASpur=[]
	arrTargetIndex=[]
	for i in arrRange_chunk:
		#For each chunk,
		## if it has "spurious exon", spurious exon and insertion should be relocated
		arrIndex=range(i[0], i[1]+1)		

		#Find AA locate within megachunk
		bIsChunkWithSpurExon=False
		listSus=list( set(arrIdx_sus).intersection(set(arrIndex)) );
		if len(listSus)>0:
			bIsChunkWithSpurExon=True;

		if bIsChunkWithSpurExon==False:
			continue

		arrQAA_inchunk=list( set( arrIdx_query).intersection(set(arrIndex)) )
		arrSAA_inchunk=list( set( arrIdx_subject).intersection(set(arrIndex) )) 
	

		strQueryExon=strSubjectExon="";
		for j in arrQAA_inchunk:
			strQueryExon=strQueryExon+a_strQuery[j]
			arrStrQ[j]="-"
			arrStrA[j]="-"
		for j in arrSAA_inchunk:
			strSubjectExon=strSubjectExon+a_strSubject[j];
			arrStrS[j]="-"
			arrStrA[j]="-"

	
		nLen_querynew=len(strQueryExon)
		nLen_subjnew=len(strSubjectExon)
		nLen_totalnew=nLen_querynew+nLen_subjnew


		nTargetIndex=i[1]	
		arrTargetIndex.append( nTargetIndex )
		arrRange_QSpur.append( "-"*(nLen_totalnew-nLen_querynew)+strQueryExon )	
		arrRange_SSpur.append( strSubjectExon+"-"*(nLen_totalnew-nLen_subjnew) )
		arrRange_ASpur.append( "-"*nLen_totalnew )

	strQuery_fix=insertbyindex( "".join(arrStrQ), arrTargetIndex, arrRange_QSpur )
	strAlign_fix=insertbyindex( "".join(arrStrA), arrTargetIndex, arrRange_ASpur)
	strSubject_fix=insertbyindex( "".join(arrStrS), arrTargetIndex, arrRange_SSpur )
	return (strQuery_fix, strAlign_fix, strSubject_fix)	


def FindAnchor( a_strAlignPerBase, a_strAlignPerBase_filt):
	#the first and the last index that are not part of trimmed region
	nLen_before=len(a_strAlignPerBase)
	nLen_after=len(a_strAlignPerBase_filt)
	arrIndexMatch=[]
	for i in range( 0, min(nLen_before, nLen_after) ):
		if a_strAlignPerBase[i]!=a_strAlignPerBase_filt[i]:		
			continue;
		if a_strAlignPerBase_filt[i]!=" ":
			arrIndexMatch.append( i );
	if len(arrIndexMatch)==0:
		return None;
	else:	
		return (min(arrIndexMatch), max(arrIndexMatch))




def MatchByCDSGcoord( a_strqueryID, a_strrefID ):
	#Using CDS information, find aminoacid coordinate where they can't be align to each other
	strStrand=dictQueryGTF[ a_strqueryID ][0]
	arrGcoord_query=dictQueryGTF[ a_strqueryID ][1]
	arrGcoord_ref=dictRefGTF[ a_strrefID ][1]


	arrGcoord_query_match=[]
	for i in arrGcoord_query:
		if i in arrGcoord_ref:
			arrGcoord_query_match.append( True )
		else:
			arrGcoord_query_match.append( False )

	arrGcoord_ref_match=[]
	for i in arrGcoord_ref:
		if i in arrGcoord_query:
			arrGcoord_ref_match.append( True )
		else:
			arrGcoord_ref_match.append( False );

	#Convert it to aminoacid array
	arrAA_query_match=[];
	nLenProt_query=int( len( arrGcoord_query_match )/3 )
	for i in range( 0, nLenProt_query ):
		nBaseMatch=sum( [arrGcoord_query_match[i*3], arrGcoord_query_match[i*3+1], arrGcoord_query_match[i*3+2]] )
		bCodonMatch=True if nBaseMatch==3 else False;
		arrAA_query_match.append( bCodonMatch );
	
	arrAA_ref_match=[]
	nLenProt_ref=int( len( arrGcoord_ref_match )/3 )
	for i in range( 0, nLenProt_ref ):
		nBaseMatch=sum( [arrGcoord_ref_match[i*3], arrGcoord_ref_match[i*3+1], arrGcoord_ref_match[i*3+2]] );
		bCodonMatch=True if nBaseMatch==3 else False;
		arrAA_ref_match.append( bCodonMatch );

	if strStrand=="-":
		arrAA_query_match.reverse();
		arrAA_ref_match.reverse();


	bHasImpossibleMatch=(sum(arrAA_query_match)!=len(arrAA_query_match)) | (sum(arrAA_ref_match)!=len(arrAA_ref_match))

	return([bHasImpossibleMatch, arrAA_query_match, arrAA_ref_match] )



def MarkWrongAlign( a_strQueryTX, a_strRefTX, a_strQuery, a_strAlignPerBase, a_strRef):		
	#0) Remove AA alignment without gcoord alignment
	#From the alignment text, what is the AA position of each qury/ref protein sequence, 1-based

	arrPos_query=[]; nPos_query=0
	arrQuery=list( a_strQuery )
	arrRef=list( a_strRef )
	for i in range(0, len(a_strQuery)):
		if a_strQuery[i].isalpha():
			nPos_query+=1;
			arrPos_query.append( nPos_query )	
		else:
			arrPos_query.append( -1 )

	arrPos_ref=[]; nPos_ref=0
	for i in range(0, len(a_strRef) ):
		if a_strRef[i].isalpha():
			nPos_ref+=1;
			arrPos_ref.append( nPos_ref )
		else:	
			arrPos_ref.append( -1 );

	arrAlignPerBase=list( a_strAlignPerBase )
	for i in range(0, len(a_strAlignPerBase)):
		if a_strAlignPerBase[ i ] not in ["|", "."]:
			continue;

		nAApos_query=arrPos_query[ i ]
		nAApos_ref=arrPos_ref[ i ]

		if nAApos_query < 0 or nAApos_ref < 0:
			print("impossible")
			continue;

		#For each alignmed AA position, find whether their CDS gcoord is overlapping or not
		arrgcoord_query=GetArrGcoordByAApos( a_strQueryTX, dictQueryGTF, nAApos_query )
		arrgcoord_ref=GetArrGcoordByAApos( a_strRefTX, dictRefGTF, nAApos_ref)
		nOverlap=len( set( arrgcoord_query ).intersection(set(arrgcoord_ref)) );

		if nOverlap==0:
			arrAlignPerBase[ i ]="&";

	strAlignPerBase_filt="".join( arrAlignPerBase )
	strQuery_filt=a_strQuery
	strRef_filt=a_strRef
	if bTest==True:
		print("\n".join(["1", strQuery_filt, strAlignPerBase_filt, strRef_filt]))

	#2) Remove alignment with too much ".". Likely chimeric ends
	nLeftSpaceCount=strAlignPerBase_filt.count(" ");
	arrChunk=strAlignPerBase_filt.lstrip(" ").split("-")	
	arrMismatchExon=[]
	for i in arrChunk:
		nLen=len(i)
		nMismatch=i.count(".")	

		if nLen==0:
			bWrongExon=False;
		else:
			bWrongExon=True if float(nMismatch)/float(nLen) >= 0.3 else False;
		arrMismatchExon.append( bWrongExon );

	if sum(arrMismatchExon)>0:		
		for i in range( 0, len(arrChunk ) ):
			if arrMismatchExon[i]==True:	
				arrChunk[ i ]="&"*len( arrChunk[ i ] )
	
		strAlignPerBase_filt=" "*nLeftSpaceCount+"-".join(arrChunk)		
	if bTest==True:
		print("\n".join(["2", strQuery_filt, strAlignPerBase_filt, strRef_filt]))

	#3) Remove spurious alignment at boundaries of alignment
	nLeftSpaceCount=strAlignPerBase_filt.count(" ");
	arrChunk=strAlignPerBase_filt.lstrip(" ").split("-")
	nWnd=10
	arrSpurBorderExon=[]
	for i in arrChunk:
		nLen=len(i)
		if nLen==0:
			bWrongExon=False;
		else:
			bWrongExon=True if i[:nWnd].find(".")>=0 or i[-nWnd:].find(".")>=0 else False;	
		arrSpurBorderExon.append(bWrongExon) 

	if sum(arrSpurBorderExon)>0:
		for i in range( 0, len(arrChunk ) ):
			if arrSpurBorderExon[i]==True:	
				strAlignChunk=arrChunk[ i ];
				nLen=len( strAlignChunk )

				nTrimLeft=nTrimRight=0			
	
				if nLen >= 2*nWnd:	
					nIdx_trimleft=nWnd-strAlignChunk[:nWnd][::-1].find(".")-1 if strAlignChunk[:nWnd].find(".") >=0 else -1;
					nIdx_trimright=strAlignChunk[-nWnd:].find(".") if strAlignChunk[-nWnd:].find(".") >=0 else nWnd;
					nTrimLeft=nIdx_trimleft+1 if nIdx_trimleft>=0 else 0
					nTrimRight=nWnd-nIdx_trimright
				else:
					nWnd_left=int(nLen/2)
					nWnd_right=nLen-nWnd_left;
					
					nIdx_trimleft=nWnd_left-strAlignChunk[:nWnd_left][::-1].find(".")-1 if strAlignChunk[:nWnd_left].find(".") >=0 else -1;
					nIdx_trimright=strAlignChunk[-nWnd_right:].find(".") if strAlignChunk[-nWnd_right:].find(".") >=0 else nWnd_right;
					nTrimLeft=nIdx_trimleft+1 if nIdx_trimleft>=0 else 0
					nTrimRight=nWnd_right-nIdx_trimright

				nIdx_left_retain=nTrimLeft;
				nIdx_right_retain=nLen-nTrimRight
			
				arrChunk[ i ]="&"*nTrimLeft+arrChunk[ i ][ nIdx_left_retain:nIdx_right_retain ]+"&"*nTrimRight
		strAlignPerBase_filt=" "*nLeftSpaceCount+"-".join(arrChunk)
	if bTest==True:
		print("\n".join(["3", strQuery_filt, strAlignPerBase_filt, strRef_filt]))
		
	#if there is no ".", nTrimRight==0

	


	#4) For all & location in gap, check the possibility of realignment
	if strAlignPerBase_filt.count("&")>0 and strAlignPerBase_filt.count("-")>0:	
		arrAlignPerBase=list( strAlignPerBase_filt )
		arrGap=strAlignPerBase_filt.split("|");
		arrCandGap=[];
		arrAllCandPos=[ i+1  for i in range(0, len(arrAlignPerBase)) if arrAlignPerBase[i]=="&"  ];    #all position of &

		bInnerAlignUpdate=False;
		for i in arrGap:
			nLen=len( i );
			if nLen==0:
				bCand=False;
			else:
				bCand=True if i.count( "&" )>0 and i.count("-") >0 else False; 
			arrCandGap.append( bCand )

		if sum( arrCandGap )>0:
			arrPosPerCandGap=[];	#position of all letters
			arrCandPosPerCand=[];	#position of & letters	
			
			nPosStart=1
			#acroos all gaps, if gap is candidate for realignment, try realignment

				
			for i in range(0, len( arrGap )):
				nLen=len( arrGap[i] )	
				if nLen==0:
					nLen=1

				nPosEnd=nPosStart+nLen-1;
				if False==arrCandGap[ i ]:
					nPosStart=nPosEnd+1;
					continue;

				if arrGap[ i ].count(" " )>0:	
					# the following procedure is only for the gap between matched alignment
					# for instance "------" at the end is not needed to undergo this step
					nPosStart=nPosEnd+1;
					continue;

				#obtain all cand position within the gap
				arrCandPosInGap=[j for j in arrAllCandPos if j>=nPosStart and j<=nPosEnd]				


				#some aminoacid alignment is very short with gap
				#it is possible where the alignment is actually better located at the other side of the gap

				#left tilting
				arrCandResolvedLeft=[];
				nRelPos=1

				nPosStart_targetsearch=nPosStart
				for j in arrCandPosInGap:
					nSourcePos=j
					#At the source position, does it have aligned AA in query/sample?
					#	    *
					#Q----------Q
					#-----------&
					#-----------Q
				
					#bGapInQuery: at the starred position, 
					
					
					for k in range( nPosStart_targetsearch, nPosEnd+1):
						nTargetPos=k;
						if nSourcePos==nTargetPos:
							continue;

						#at the target position, which of query/ref has the gap?
						bGapInQuery=True if a_strQuery[ nTargetPos-1 ]=="-" else False;
						nAApos_from=arrPos_query[ nSourcePos-1 ] if bGapInQuery else arrPos_ref[ nSourcePos-1 ]
						nAApos_to=arrPos_ref[ nTargetPos-1 ] if bGapInQuery else arrPos_query[ nTargetPos-1 ]
						
						arrgcoord_from=GetArrGcoordByAApos( a_strQueryTX, dictQueryGTF, nAApos_from ) if bGapInQuery else GetArrGcoordByAApos( a_strRefTX, dictRefGTF, nAApos_from )
						arrgcoord_to=GetArrGcoordByAApos( a_strRefTX, dictRefGTF, nAApos_to ) if bGapInQuery else GetArrGcoordByAApos( a_strQueryTX, dictQueryGTF, nAApos_to )

						chrFrom=a_strQuery[ nSourcePos-1 ] if bGapInQuery else a_strRef[ nSourcePos-1 ]
						chrTo=a_strRef[ nTargetPos-1 ] if bGapInQuery else a_strQuery[ nTargetPos-1 ];
						if chrFrom!=chrTo:
							break;		
						nLenMatch=len( set( arrgcoord_from ).intersection( set(arrgcoord_to)) )
						if nLenMatch==0:
							break;

						nPosStart_targetsearch+=1;
						#if "match" status at the particular "source" location should be moved to the "target" location,
						#we need to update the query sequence and ref sequence
						arrAlignPerBase[ nSourcePos-1 ]="U" if bGapInQuery else "D"	#at source location, delete up (query), delete down (ref)
						arrAlignPerBase[ nTargetPos-1 ]="Q" if bGapInQuery else "R"	#at target location, add query, add ref
						arrCandResolvedLeft.append( j )

				arrPosCandPos=list( set( arrCandPosInGap )-set( arrCandResolvedLeft )	)
	
				##right tilting	
				arrCandResolvedRight=[];
				nRelPos=len( arrPosCandPos )
				nPosEnd_targetsearch=nPosEnd
				for j in arrPosCandPos[::-1]:
					bIsOver=False;
					for k in range( nPosStart, nPosEnd_targetsearch+1)[::-1]:	
						nTargetPos=k;
						if nSourcePos==nTargetPos:
							continue;

						#at the target position, which of query/ref has the gap?
						bGapInQuery=True if a_strQuery[ nTargetPos-1 ]=="-" else False;

						nAApos_from=arrPos_query[ nSourcePos-1 ] if bGapInQuery else arrPos_ref[ nSourcePos-1 ]
						nAApos_to=arrPos_ref[ nTargetPos-1 ] if bGapInQuery else arrPos_query[ nTargetPos-1 ]
						
						arrgcoord_from=GetArrGcoordByAApos( a_strQueryTX, dictQueryGTF, nAApos_from ) if bGapInQuery else GetArrGcoordByAApos( a_strRefTX, dictRefGTF, nAApos_from )
						arrgcoord_to=GetArrGcoordByAApos( a_strRefTX, dictRefGTF, nAApos_to ) if bGapInQuery else GetArrGcoordByAApos( a_strQueryTX, dictQueryGTF, nAApos_to )
			
						#character
						chrFrom=a_strQuery[ nSourcePos-1 ] if bGapInQuery else a_strRef[ nSourcePos-1 ]
						chrTo=a_strRef[ nTargetPos-1 ] if bGapInQuery else a_strQuery[ nTargetPos-1 ];
						if chrFrom!=chrTo:
							bIsOver=True;
							break;

						nLenMatch=len( set( arrgcoord_from ).intersection( set(arrgcoord_to)) )
						if nLenMatch==0:
							bIsOver=True
							break;
										
						nPosEnd_targetsearch-=1
						arrAlignPerBase[ nSourcePos-1 ]="U" if bGapInQuery else "D"     #at source location, delete up (query), delete down (ref)
						arrAlignPerBase[ nTargetPos-1 ]="Q" if bGapInQuery else "R"     #at target location, add query, add ref
						arrCandResolvedLeft.append( j )					
					if bIsOver:
						break;		
				arrPosCandPos=list( set( arrPosCandPos )-set( arrCandResolvedRight ) )		
				nPosStart=nPosEnd+1;

				if len(arrCandResolvedLeft)>0 or len(arrCandResolvedRight) > 0:
					bInnerAlignUpdate=True;

		strAlignPerBase_filt="".join( arrAlignPerBase )

		
		if bInnerAlignUpdate:			
			if bTest==True:
				print(strAlignPerBase_filt)

			#update all alignmet track
			arrQuery=list( strQuery_filt )
			arrRef=list( strRef_filt )
			for i in range( 0, len(strAlignPerBase_filt) ):
				if strAlignPerBase_filt[i]=="U":
					arrQuery[i]="-"	
				elif strAlignPerBase_filt[i]=="D":
					arrRef[i]="-"
				elif strAlignPerBase_filt[i]=="Q":
					arrQuery[i]=arrRef[i]
				elif strAlignPerBase_filt[i]=="R":
					arrRef[i]=arrQuery[i]
				else:
					continue;
			strQuery_filt="".join(arrQuery)
			strRef_filt="".join(arrRef)
			strAlignPerBase_filt=strAlignPerBase_filt.replace("U", "-").replace("D", "-")
			strAlignPerBase_filt=strAlignPerBase_filt.replace("Q", "|").replace("R", "|")

	if bTest==True:
		print("\n".join(["4", strQuery_filt, strAlignPerBase_filt, strRef_filt]))		


	#5) Replace "^, &, *" to "-"
	strAlignPerBase_filt=strAlignPerBase_filt.replace("&", "-")
	strAlignPerBase_filt=strAlignPerBase_filt.rstrip("-")	
	

	nLeftSpaceCount=strAlignPerBase_filt.count(" ");
	strAlignPerBase_filt_noleftspace=strAlignPerBase_filt.lstrip(" ");
	nLeftStripDash=len( strAlignPerBase_filt_noleftspace )-len( strAlignPerBase_filt_noleftspace.lstrip("-" ));
	strAlignPerBase_filt=" "*(nLeftSpaceCount+nLeftStripDash)+strAlignPerBase_filt_noleftspace.lstrip("-")

	if bTest==True:
		print("\n".join(["5", strQuery_filt, strAlignPerBase_filt, strRef_filt]))

	#6) Remove dash
	arrResult=RemoveDash( strQuery_filt, strAlignPerBase_filt, strRef_filt)
	if bTest==True:
		print("\n".join(["6"]+arrResult ));		
	return( arrResult )




	
def GetStringNoDash(a_strInput, a_start, a_end ):
	strInputLeft=a_strInput[:a_start]
	strInputMatch=a_strInput[a_start:a_end]
	strInputRight=a_strInput[a_end:]
	
	nDashLeft=strInputLeft.count("-")
	nDashRight=strInputRight.count("-")
	
	if bTest:
		print("GetStringNoDash")
		print([a_strInput, a_start, a_end])
		print([strInputLeft, strInputMatch, strInputRight])
		
		
	return( " "*nDashLeft+"".join( strInputLeft.split("-") )+strInputMatch+"".join( strInputRight.split("-") ) );

	




def RemoveDash( a_strQuery, a_strUpdateAlignBase, a_strRef ):
	##As a result of removal of spurious alignment, some query/ref sequence may have "-" in query/ref without align information
	##
##	return( [a_strQuery, a_strRef] )

	#Using the updatealignbase, I can find the first position of "alignment" denoted as "|"
	## nIndex_nonblank: the first position of "|"
	
	nIndex_nonblank=len(a_strUpdateAlignBase)-len( a_strUpdateAlignBase.lstrip(" ") )
	nIndex_nonblank=0 if nIndex_nonblank<0 else nIndex_nonblank


	strQueryUpdate=GetStringNoDash( a_strQuery, nIndex_nonblank, len(a_strUpdateAlignBase))
	strRefUpdate=GetStringNoDash( a_strRef, nIndex_nonblank, len(a_strUpdateAlignBase))
	

	#Remove unused space
	nSpace_query=len( strQueryUpdate )-len(strQueryUpdate.lstrip(" "));
	nSpace_ref=len( strRefUpdate )-len(strRefUpdate.lstrip(" "))
	nSpace_alignbase=len( a_strUpdateAlignBase )-len( a_strUpdateAlignBase.lstrip(" "));
	nSpace_min=min( nSpace_query, nSpace_ref, nSpace_alignbase )


	return( [strQueryUpdate[nSpace_min:], a_strUpdateAlignBase[nSpace_min:], strRefUpdate[nSpace_min:] ] )
 







print("Fix spurious alignment");
nIndex=0
strHeader=""
strQueryID=strRefID=""
for line in ifalign:
	strLine=line.rstrip("\n");
	if (nIndex % 5)==0:
		arr=strLine.split(",")
		strQueryID=arr[0]	#protein name
		strQueryTXs=arr[1]	#list of transcripts
		strRefID=arr[2]		#best matching ref transcript id
		strHeader=strLine	

		strQueryRepTX=strQueryTXs.split(";")[0]
		strQueryRepTX=strQueryRepTX.split("(")[0]
	elif (nIndex % 5)==1:
		strQuery=strLine
	elif (nIndex % 5)==2:
		strAlignPerBase=strLine;
	elif (nIndex % 5)==3:
		strSubject=strLine
	else:
		#if strQueryID=="PB.22532.3":
		#print(strQueryID)
		#print("\n".join([strQuery, strAlignPerBase, strSubject])+"\n")
		#print("-------")
	#	if strQueryID !="ENSG00000147099.19_prot_2":
	#		nIndex=nIndex+1;
		#	continue;
		if bTest==True:
			if strQueryID!="ENSG00000100075.9_prot_2":
				nIndex=nIndex+1;
				continue;

		if strRefID=="NA":
			strQuery_result=strAlignPerBase_result=strSubject_result="empty"
		else:
			#0) Peptides from difference exons may align each other. which is false positive alignment	
			## Based on CDS comparison, some peptide of query/ref may be unique to either query/ref
			## In some cases, they align each other with many mismatches and gap (false positive)
			## step 0 is to eliminate those error prone alignment

		
			arrFinal=MarkWrongAlign( strQueryRepTX, strRefID, strQuery, strAlignPerBase, strSubject )


			strQuery_result=arrFinal[ 0 ]
			strSubject_result=arrFinal[ 2 ]
			strAlignPerBase_result=arrFinal[1 ]
			#	strAlignPerBase_result="empty"

		if bTest==True:
			print(strHeader)
			print("\n".join( [strQuery_result, strAlignPerBase_result, strSubject_result] )+"\n\n")	

		ffile.write(strHeader+"\n")
		ffile.write("\n".join( [strQuery_result, strAlignPerBase_result, strSubject_result] )+"\n\n");	

	nIndex=nIndex+1;

ffile.close()
	








	
	




