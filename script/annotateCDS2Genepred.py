import sys
import copy
import pysam

igenepred=sys.argv[1]
iCDS=sys.argv[2]
ogenepred=sys.argv[3]

#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd
#TALONT000093124	chr1	-	184874	199900	185219	187832	

#exonCount	exonStarts	exonEnds	score
#11	184874,185490,186316,187128,187375,187754,188129,188438,188790,195262,199836,	185350,185559,186469,187287,187577,187890,188266,188584,188902,195416,199900,	0

#genename	cdsStartStat	cdsEndStat exonFrames
#TALONT000093124	none	none	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,


dictCDS_coord={}
with pysam.FastxFile( iCDS ) as fin:
	for entry in fin:
		arr=entry.name.split("|")
		strTxname=arr[0].split(":")[0];
		
		arr2=entry.comment.split(" ")
		arrCDScoord=arr2[3].split(":")[1].split("-")


		#ANGEL output has 1-based coordinate
		##ANGEL's nEnd includes stop codon
		##Genepred CDS stop also incluids stop codon
		##However, GTF does not
		nStart=int( arrCDScoord[0] )
		nEnd=int( arrCDScoord[1] )
		dictCDS_coord[ strTxname ]=(nStart, nEnd)

print("CDS")
print("loading CDS gcoord")

def GetArrExon( a_exonstart, a_exonend ):
	#arrStart, arrEnd is 1-based coordinate
	arrStart=[int(i)+1 for i in a_exonstart.rstrip(",").split(",")];
	arrEnd=[int(i) for i in a_exonend.rstrip(",").split(",")];
		
	arrExon=[]
	for i in range(0, len(arrStart)):
		arrExon=arrExon+[(arrStart[i], arrEnd[i])]
	return( arrExon )

def GetArrCDS( a_arrexon, a_cdsstart, a_cdsend):
	arrCDS=[]
	for i in a_arrexon:
		if i[1]<a_cdsstart:
			continue;
		if i[0]>a_cdsend:
			continue;

		nCDSstart=max( a_cdsstart, i[0])
		nCDSend=min( a_cdsend, i[1]);
	
		arrCDS=arrCDS+[(nCDSstart, nCDSend)]	
		
	return(arrCDS)


def GetArrGcoordExon(  a_arrexon ):
	#arrStart, arrEnd is 1-based coordinate
	arrGcoord=[]    #1-based coordinate of all exons
	arrIndexStart=[];       #index of start in arrGcoord
	for i in a_arrexon :
		arrIndexStart.append( i[0] );
		arrGcoord=arrGcoord+list(range( i[0], i[1]+1 ));
	return (arrGcoord, arrIndexStart);

def GetArrGcoordCDS( a_arrGcoord, a_idxstart, a_idxend, a_strand):	
	arrGcoordCDS=[];
	if a_strand=="+":
		arrGcoordCDS=a_arrGcoord[ (a_idxstart-1):(a_idxend) ]
	else:
		arrGcoordCDS=a_arrGcoord[ (-a_idxend):]
		arrGcoordCDS=arrGcoordCDS[:(a_idxend-a_idxstart+1)]
	return( arrGcoordCDS )


def GetORF( a_mockorf, a_tpCDScoord, a_arrIndexStart, a_nExonLen, a_strand ):
	nBase_tillStart=a_tpCDScoord[0]-1 if a_strand=="+" else a_nExonLen-a_tpCDScoord[1]        #number of bases before start codon	
	nNameoji=(nBase_tillStart % 3)
	nOffset=3-nNameoji if a_strand=="+" else 3+nNameoji;
	arrGcoordORF=a_mockorf[ nOffset:(nOffset+a_nExonLen) ]
	arrORF=[ str(a_mockorf[i]) for i in a_arrIndexStart ]
	return(arrORF)


def GetCDSCoord( a_exonstart, a_exonend, a_tpCDScoord, a_strand ):	
	#a_tpCDScoord: 1-based coordinate of CDS against mRNA
	#		[0]: start coordinate of CDS
	#		[1]: end coordinate of stop codon which is outside of CDS

	nIdxStart=a_tpCDScoord[0]
	nIdxEnd=a_tpCDScoord[1]
	nCDSLen=nIdxEnd-nIdxStart+1	
	arrExon=GetArrExon( a_exonstart, a_exonend)


	#1-based coordinate
	tpArrGcoordExon=GetArrGcoordExon( arrExon );
	arrGcoord=tpArrGcoordExon[0]
	arrIndexStart=tpArrGcoordExon[1]
	nExonLen=len(arrGcoord)
	
	#1-based coordinate
	arrGcoordCDS=GetArrGcoordCDS( arrGcoord, nIdxStart, nIdxEnd, a_strand )
	arrCDS=GetArrCDS( arrExon, arrGcoordCDS[0], arrGcoordCDS[1])

	arrORF=[]
	exonNumbers=list( range(0, len(arrExon)));
	cdsNumbers=list( range(0, len(arrCDS)) );
	if a_strand=="-":
		exonNumbers.reverse();
		cdsNumbers.reverse();
	
	nCumCDS=0;
	for i in exonNumbers:
		setBaseInExon=set( range( arrExon[i][0], arrExon[i][1]+1 ) );
		
		#if exon has no CDS
		nCDSsize=len( setBaseInExon.intersection( set(arrGcoordCDS) ) )	
		if nCDSsize==0:
			arrORF=arrORF+[-1];
			continue;
		
		#if exon has CDS
		if i==0 and a_strand=="+":
			arrORF=arrORF+[0];
		elif i==(max( cdsNumbers )-1) and a_strand=="-":
			arrORF=arrORF+[0];
		else:
			arrORF=arrORF+[nCumCDS % 3];
			
		nCumCDS=nCumCDS+nCDSsize;

	if a_strand=="-":
		arrORF.reverse();
	arrORF=[str(i) for i in arrORF];		
	return (arrGcoordCDS[0], arrGcoordCDS[-1], arrORF);	

#Find Gcoord using genepred input
ofile=open(ogenepred, "w+");
for line in open(igenepred):
	arr=line.rstrip("\n").split("\t");

	txid=arr[0]
	if txid in dictCDS_coord:
#		if txid not in ["TALONT000093124"]:# ["TALONT000095358", "TALONT000096303", "TALONT000102070"]:
#			continue;

		#print(txid)
		tpCDS=GetCDSCoord( arr[8], arr[9], dictCDS_coord[txid], arr[2])
	
		arr[5]=str( tpCDS[0]-1 )
		arr[6]=str( tpCDS[1] )

		arrORF=tpCDS[2]+[""]
		arr[14]=",".join(arrORF);		

#		if txid== "TALONT000567675":
#			print(tpCDS)
		#print("\t".join(arr))	
		ofile.write("\t".join(arr)+"\n");
	else:
		ofile.write("\t".join(arr)+"\n");
ofile.close();







