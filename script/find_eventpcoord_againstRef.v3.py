import sys
import re

ifpath=sys.argv[1]
ofpath=sys.argv[2]

bIgnoreINDELMISMATCH=False;

##Convert the pcoort of event against PB protein to Ref protein
#ifpath="/scratch/jmaeng/longreadCAGE/h1299/analysis_hg38_noalt.v2/ORFprediction/wGTF.v2/LRCAGE.long_ORF.framecall_align.txt"
#ifpath="/scratch/jmaeng/longreadCAGE/h1299/analysis_hg38_noalt.v2/peptidome_analysis.v2/alternative_splicing/test_input.txt"
#ofpath="/scratch/jmaeng/longreadCAGE/h1299/analysis_hg38_noalt.v2/peptidome_analysis.v2/alternative_splicing/ref"
def ConvertList2Tuple(a_arr):
	arrGroup=[]
	if len( (a_arr) )>0:
		arrIdx_notconsecutive=[ i for i in range(0, len( a_arr )-1) if a_arr[i+1]-a_arr[i]>1] + [len(a_arr)-1];
		for i in range( 0, len(arrIdx_notconsecutive) ):
			if i==0:
				arrGroup.append( [ a_arr[ 0 ], a_arr[ arrIdx_notconsecutive[i] ] ]);
			else:	
				arrGroup.append( [ a_arr[ arrIdx_notconsecutive[i-1]+1 ], a_arr[ arrIdx_notconsecutive[i] ] ]);
	
	return( arrGroup )


def GetArrPcoord( a_str):
	#alignbase to pcoord (1-based)
	nIndex=0;
	arr=[];
	for i in a_str:
		if i.isalpha():
			nIndex=nIndex+1;
			arr.append( nIndex)
		else:
			#0 is place holder for non alphabet base
			arr.append(0)
	return arr


def GetVARgroup( a_strQuery, a_strAlignPerBase, a_strSubject ):
	#0-based coordinate
	arrIdx_gap_inquery=[x.start() for x in re.finditer("-", a_strQuery)]
	arrIdx_gap_inalign=[x.start() for x in re.finditer("-", a_strAlignPerBase)]
	arrIdx_gap_insubject=[x.start() for x in re.finditer("-", a_strSubject)]	
	
	arrIdx_nonempty_query=[x.start() for x in re.finditer("[A-z,-]", a_strQuery)]
	arrIdx_nonempty_align=[x.start() for x in re.finditer("[|,.,-]", a_strAlignPerBase)]
	arrIdx_nonempty_subject=[x.start() for x in re.finditer("[A-z,-]", a_strSubject)]

	arrPCoordPair=[]; #pcoord against query, pcoord against ref
	if len(arrIdx_nonempty_align)==0:
		arrPCoordPair.append( [ "OUTFRAME",[-1, -1, 0], [-1, -1, 0]] );
		return arrPCoordPair

	nIndex_nonempty_match_first=min(arrIdx_nonempty_align)
	nIndex_nonempty_match_last=max(arrIdx_nonempty_align)

	arrIdx_chimer=list( set(arrIdx_nonempty_query)-set(arrIdx_nonempty_align) ); arrIdx_chimer.sort();
	arrIdx_chimer_ref=list( set(arrIdx_nonempty_subject)-set(arrIdx_nonempty_align) ); arrIdx_chimer_ref.sort();

	arrIdx_chimer5=[i for i in arrIdx_chimer if i<nIndex_nonempty_match_first]
	arrIdx_chimer3=[i for i in arrIdx_chimer if i>nIndex_nonempty_match_last]	
	arrIdx_chimer_ref5=[i for i in arrIdx_chimer_ref if i<nIndex_nonempty_match_first]
	arrIdx_chimer_ref3=[i for i in arrIdx_chimer_ref if i>nIndex_nonempty_match_last]


	#pcoord of query, sample protein
	#0-based coordinate
	arrIdx_chimer5_group=ConvertList2Tuple( arrIdx_chimer5 )
	arrIdx_chimer3_group=ConvertList2Tuple( arrIdx_chimer3 )
	arrIdx_chimer_ref5_group=ConvertList2Tuple( arrIdx_chimer_ref5)
	arrIdx_chimer_ref3_group=ConvertList2Tuple( arrIdx_chimer_ref3)


	#if you put the index of each alignperbase coord, you can get pcoord
	arrToPcoord_query=GetArrPcoord(a_strQuery)
	arrToPcoord_subject=GetArrPcoord(a_strSubject)


	#convert to AA coord
	for i in arrIdx_chimer5_group:
		nIndex_start=i[0]; nIndex_end=i[1]
		nstart_qpcoord=arrToPcoord_query[ nIndex_start ]
		nend_qpcoord=arrToPcoord_query[ nIndex_end ]
		nstart_spcoord=arrToPcoord_subject[ nIndex_end+1 ]
		nend_spcoord=arrToPcoord_subject[ nIndex_end+2 ]			
		arrPCoordPair.append( [ "CHIMERIC_5", [nstart_qpcoord, nend_qpcoord, nend_qpcoord-nstart_qpcoord+1], [nstart_spcoord, nend_spcoord, nend_spcoord-nstart_spcoord+1] ] )	

	for i in arrIdx_chimer3_group:
		nIndex_start=i[0]; nIndex_end=i[1]
		nstart_qpcoord=arrToPcoord_query[ nIndex_start ]
		nend_qpcoord=arrToPcoord_query[ nIndex_end ]
		nstart_spcoord=arrToPcoord_subject[ nIndex_start-2 ]
		nend_spcoord=arrToPcoord_subject[ nIndex_start-1 ]
		arrPCoordPair.append( [ "CHIMERIC_3", [nstart_qpcoord, nend_qpcoord, nend_qpcoord-nstart_qpcoord+1], [nstart_spcoord, nend_spcoord, nend_spcoord-nstart_spcoord+1] ] )

	for i in arrIdx_chimer_ref5_group:
		nIndex_start=i[0]; nIndex_end=i[1]
		nstart_qpcoord=arrToPcoord_query[ nIndex_end+1 ]
		nend_qpcoord=arrToPcoord_query[ nIndex_end+1 ]
		nstart_spcoord=arrToPcoord_subject[ nIndex_start ]
		nend_spcoord=arrToPcoord_subject[ nIndex_end]
		arrPCoordPair.append( [ "TRUNCATED_5", [nstart_qpcoord, nend_qpcoord, nend_qpcoord-nstart_qpcoord+1], [nstart_spcoord, nend_spcoord, nend_spcoord-nstart_spcoord+1] ] )

	for i in arrIdx_chimer_ref3_group:
		nIndex_start=i[0]; nIndex_end=i[1]
		nstart_qpcoord=arrToPcoord_query[ nIndex_start-1 ]	
		nend_qpcoord=arrToPcoord_query[ nIndex_start-1 ]
		nstart_spcoord=arrToPcoord_subject[ nIndex_start ]
		nend_spcoord=arrToPcoord_subject[ nIndex_end]
		arrPCoordPair.append( [ "TRUNCATED_3", [nstart_qpcoord, nend_qpcoord, nend_qpcoord-nstart_qpcoord+1], [nstart_spcoord, nend_spcoord, nend_spcoord-nstart_spcoord+1] ] )

	if False==bIgnoreINDELMISMATCH:
		arrIdx_del=list( set( arrIdx_gap_inquery ) & set( arrIdx_gap_inalign ) ); arrIdx_del.sort();
		arrIdx_ins=list( set( arrIdx_gap_inalign ) - set( arrIdx_gap_inquery ) ); arrIdx_ins.sort();
		arrIdx_mismatch=[x.start() for x in re.finditer("[.]", a_strAlignPerBase)]
		arrIdx_ins_group=ConvertList2Tuple( arrIdx_ins )
		arrIdx_del_group=ConvertList2Tuple( arrIdx_del )

		for i in arrIdx_ins_group:
			#   ***********         coordinate against PB
			#BBBAAAAAAAAAAABBBB
			#|||-----------||||
			#BBB-----------BBBB
			#  *           *        coordinate against Ref
			
			nIndex_start=i[0]; nIndex_end=i[1]
			#nstart_qpcoord, nstart_spcoord 1-based coordinate
			nstart_qpcoord=arrToPcoord_query[ nIndex_start ]
			nend_qpcoord=arrToPcoord_query[ nIndex_end ]+1
	
			#Insertion may provide R, K which provides tryptic position. Therefore, I extend end cooridnate by 1
			nstart_spcoord=arrToPcoord_subject[ nIndex_start-1 ]
			nend_spcoord=arrToPcoord_subject[ nIndex_end+1 ]
			
			arrPCoordPair.append( [ "INSERTION", [nstart_qpcoord, nend_qpcoord, nend_qpcoord-nstart_qpcoord+1], [nstart_spcoord, nend_spcoord, nend_spcoord-nstart_spcoord+1] ] )

		for i in arrIdx_del_group:
			#  *           *        coordinate against PB
			#BBB-----------BBBB
			#|||-----------||||
			#BBBAAAAAAAAAAABBBB
			#   ***********         coordinate against Ref
			
			nIndex_start=i[0]; nIndex_end=i[1]
			nstart_qpcoord=arrToPcoord_query[ nIndex_start-1 ]
			nend_qpcoord=arrToPcoord_query[ nIndex_end+1 ]
			
			nstart_spcoord=arrToPcoord_subject[ nIndex_start ]
			nend_spcoord=arrToPcoord_subject[ nIndex_end ]
			arrPCoordPair.append( [ "DELETION", [nstart_qpcoord, nend_qpcoord, nend_qpcoord-nstart_qpcoord+1], [nstart_spcoord, nend_spcoord, nend_spcoord-nstart_spcoord+1] ] )

		for i in arrIdx_mismatch:
			nIndex_start=i; nIndex_end=i
			nstart_qpcoord=arrToPcoord_query[ nIndex_start ]
			nend_qpcoord=arrToPcoord_query[ nIndex_end ]
			nstart_spcoord=arrToPcoord_subject[ nIndex_start ]
			nend_spcoord=arrToPcoord_subject[ nIndex_end ]
			arrPCoordPair.append( [ "MISMATCH", [nstart_qpcoord, nend_qpcoord, nend_qpcoord-nstart_qpcoord+1], [nstart_spcoord, nend_spcoord, nend_spcoord-nstart_spcoord+1] ] )

	for i in arrPCoordPair:

		if i[0] in ["CHIMERIC_5", "INSERTION", "DELETION", "MISMATCH" ]:
			#addition of novel sequence may affect the downstream sequence to be trypsinized
			i[1][1]=i[1][1]+1
			i[2][1]=i[2][1]+1

		elif i[0] == "TRUNCATED_5":
			#due to the initiator methionin truncation, the first AA is eliminated
			i[1][1]=i[1][1]+1
			i[2][1]=i[2][1]+1
		elif i[0] in ["CHIMERIC_3"]:
			#addition of novel sequence may impaact the upstream sequence to be trypsinzed prematually
			i[1][0]=i[1][0]-1;
			i[2][0]=i[2][0]-1

	return arrPCoordPair;

#
strPBID=strGENCODEID=strQuery=strSubject=strAlignPerBase=""
nIndex=0
bHasEvent=False;
ofile=open(ofpath+".diff.bed", "w+")
strType=""
nLen_pep=0
ofile_stat=open(ofpath+".stat.txt", "w+")
ofile_stat.write("\t".join(["readname", "ref_txID", "class", "class_3end", "len_read", "len_reftx", "match", "mismatch", "insertion", "deletion", "chimeric5", "chimeric3"])+"\n")

def GetClass(a_arr):
	if len(a_arr)==0:
		return ("NORMAL", "NORMAL");
	
	if a_arr[0] in ["UNANNOTATED_CDS", "OUTFRAME"]:
		return (a_arr[0], a_arr[0]);

	strEnd5="NORMAL"
	strEnd3="NORMAL";
	if "CHIMERIC_5" in a_arr and "TRUNCATED_5" in a_arr:
		strEnd5="CHIMERIC_TRUNCATED_5";
	elif "CHIMERIC_5" in a_arr:
		strEnd5="CHIMERIC_5";
	elif "TRUNCATED_5" in a_arr:
		strEnd5="TRUNCATED_5";

	if "CHIMERIC_3" in a_arr and "TRUNCATED_3" in a_arr:
		strEnd3="CHIMERIC_TRUNCATED_3";
	elif "CHIMERIC_3" in a_arr:
		strEnd3="CHIMERIC_3";
	elif "TRUNCATED_3" in a_arr:
		strEnd3="TRUNCATED_3";	

	return (strEnd5, strEnd3)


for line in open(ifpath):
	strLine=line.rstrip("\n");


	if (nIndex % 5)==0:
		#reset
		strPBID=strGENCODEID=strQuery=strSubject=strAlignPerBase=strType=""
		nLen_pep=0;
		bHasEvent=False;

		arr=strLine.split(",")
		strPBID=arr[0]; 
		strPBTX=arr[1]	
		strGENCODEID=arr[2]
		nLen_pep=int( arr[3].split(" ")[1] )
		nLen_ref=int( arr[4].split(" ")[1] )	

		#PB.10823_prot_1,PB.10823.17,ENST00000358389.7,LenQuery: 770,LenRef: 753

	elif (nIndex % 5)==1:
		strQuery=strLine
	elif (nIndex % 5)==2:
		strAlignPerBase=strLine;
	elif (nIndex % 5)==3:
		strSubject=strLine
	else: 
#		if strPBID!="ENSG00000146083.11_prot_1":	
#			nIndex=nIndex+1;
#			continue;

		nMatch=nMisMatch=nIns=nDel=nChim5=nChim3=0;
		arrClass=[]
		if strGENCODEID =="NA":
			#UNANNOTATED CDS
			strType="UNANNOTATED_CDS"
			strQuery_result="\t".join([strPBID, str(0), str(nLen_pep), strType, str(nLen_pep), "+"])
			strSubject_result="\t".join([strGENCODEID, str(-1), str(-1), strType, str(-1), "+"]);
			ofile.write(strQuery_result+"\t"+strSubject_result+"\n");	
			arrClass.append( strType );
		else:
			nMatch=strAlignPerBase.count("|")

			arrVAR=GetVARgroup( strQuery, strAlignPerBase, strSubject)	
			for i in arrVAR:
				strType=i[0]
				tuplPcoord_q=i[1]
				tuplPcoord_s=i[2]
				nLen_q=tuplPcoord_q[2]
				nLen_s=tuplPcoord_s[2]

				if strType=="OUTFRAME":
					strQuery_result="\t".join([strPBID, str(0), str(nLen_pep), strType, str(nLen_pep), "+"])
					strSubject_result="\t".join([strGENCODEID, str(-1), str(-1), strType, str(-1), "+"]);		
				else:		
					#1-based to 0-based coordinate
					strSubject_result="\t".join([strGENCODEID, str(tuplPcoord_s[0]-1), str(tuplPcoord_s[1]), strType, str(nLen_s), "+"])
					strQuery_result="\t".join([strPBID, str(tuplPcoord_q[0]-1), str(tuplPcoord_q[1]), strType, str(nLen_q), "+"])
			
					if strType=="MISMATCH":
						nMisMatch=nMisMatch+nLen_q;
					elif strType=="INSERTION":
						nIns=nIns+nLen_q
					elif strType=="DELETION":
						nDel=nDel+nLen_s;
					elif strType=="CHIMERIC_5":
						nChim5=nChim5+nLen_q	
					elif strType=="CHIMERIC_3":
						nChim3=nChim3+nLen_q

				ofile.write(strQuery_result+"\t"+strSubject_result+"\n");
				arrClass.append( strType )
		
		tuplClass=GetClass( arrClass )	
		ofile_stat.write("\t".join([strPBID, strGENCODEID, tuplClass[0], tuplClass[1], str(nLen_pep), str(nLen_ref), str(nMatch), str(nMisMatch), str(nIns), str(nDel), str(nChim5), str(nChim3)])+"\n")	
	nIndex=nIndex+1;
ofile_stat.close()
ofile.close();




