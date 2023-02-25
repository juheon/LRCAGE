
import sys
ibed=sys.argv[1]
obed=sys.argv[2]
#ibed="/scratch/jmaeng/longreadCAGE/h1299/pacbio/peak_01202020/annotsummit/tmp2.bed"
#obed="/scratch/jmaeng/longreadCAGE/h1299/pacbio/peak_01202020/annotsummit/tmp3.bed"

ifile=open(ibed, "r+");
ofile=open(obed, "w+");


class cPeak:	
	arr=[];	

	def __init__(self, a_arr):
		self.arr=a_arr;
	
	def GetArr(self):
		return( self.arr)

	def GetArr_simple(self):
		#Only output summit position, summit read count from CTSS BED6
		return( self.arr[0:(g_nColCount-7)]+[self.arr[g_colIdxSummit]]+[self.arr[g_colIdxSummitCnt ] ] );		

def GetSummit(a_arr):
	nArrLen=len(a_arr);
	if nArrLen%2==1:
		idx=int(nArrLen/2);
		return( a_arr[idx].GetArr_simple() );
	else:	#even number
		idx=int(nArrLen/2)-1;

		#input bed columns
                arr_summit=a_arr[idx].GetArr_simple();
		arr_summit2=a_arr[idx+1].GetArr_simple();

		#mid position
		nStart=int((int(arr_summit[-2])+int(arr_summit2[-2]))/2)

		arr_final=arr_summit[ 0:(g_nColCount-7) ]+[str(nStart)]+[arr_summit[-1]]

		#nStart=int( ( a_arr[ idx ].GetArr() )[ g_colIdxSummit ] )+int( ( a_arr[ idx+1 ].GetArr() )[ g_colIdxSummit ] );
		#nStart=int(nStart/2)

		#input bed columns 
		#arr_summit=a_arr[idx].GetArr_simple();
		#arr_final=(a_arr[ idx ].GetArr())[0:(g_nColCount-7)]+[str( nStart ), str( nStart+1)]+(a_arr[ idx ].GetArr())[9:13]
		return( arr_final )



#	else:	
#		for i in a_arr:
#			nPos=nPos+int( (i.GetArr()) [7] );
#		nPos=int( nPos/len( a_arr ) );
#	
#	arr_final=(a_arr[0].GetArr())[0:7]+[str(nPos), str(nPos+1)]+(a_arr[0].GetArr())[10:13]
#	return(arr_final)



nOldCntMax=-1;
arrOldPeak=[];

##If one CTSS has highest read count, report its coordinate
##If there is even number of tied, average of all tied position is returned
##            odd                  tied position in the middle is returned

##   0	1	2	3		4	5	6	7	8			9	10	11	12	13	14	15
##chr1	29344	29366	cluster4460	126	-	126	93	0.738095238095238	chr1	29344	29345	CTSS	2	-	1
##chr1	29344	29366	cluster4460	126	-	126	93	0.738095238095238	chr1	29345	29346	CTSS	1	-	1

g_nColCount=-1;
g_colIdxSummit=g_colIdxSummitCnt=-1
for line in ifile:
	arr=line.rstrip("\n").split("\t");

	if g_nColCount<0:
		#column number is read from the input bed file
		g_nColCount=len(arr);
		g_colIdxSummit=g_nColCount-6
		g_colIdxSummitCnt=g_nColCount-3

	curPeak=cPeak( arr );
	nCurCnt=int( arr[ g_colIdxSummitCnt ] );

	#initialization
	if len(arrOldPeak)==0:
		arrOldPeak=[ curPeak ];
		nOldCntMax=nCurCnt;
		continue;

	if (arrOldPeak[-1].GetArr())[0:6]==( curPeak.GetArr() )[0:6]:
		if nCurCnt<nOldCntMax:
			continue;
		elif nCurCnt==nOldCntMax:
			arrOldPeak.append( curPeak );
		else:
			arrOldPeak=[ curPeak ];
			nOldCntMax=nCurCnt;
	else:
		#Print to file
		out=GetSummit( arrOldPeak );
		#ofile.write( "\t".join( out[0:6]+[out[7]]+[out[10]] )+"\n");
		ofile.write( "\t".join( out )+"\n");
		arrOldPeak=[ curPeak ];
		nOldCntMax=nCurCnt;

out=GetSummit( arrOldPeak );
#ofile.write( "\t".join( out[0:6]+[out[7]]+[out[10]] )+"\n");
ofile.write( "\t".join( out )+"\n");
	




ifile.close();
ofile.close();

