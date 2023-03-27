import pysam
import sys


tabixfile = pysam.TabixFile(sys.argv[1]);
ofpath=sys.argv[2];


#tabixfile=pysam.TabixFile("/bar/jmaeng/genomes/hg38/gencode/gencode.v29.basic.annotation.sort.gtf.gz")
#ofpath="/scratch/jmaeng/longreadCAGE/h1299/analysis_hg38_noalt.v2/ORFprediction/NMD/gencode.v29.basic.annotation.sort.NMD.txt"

#Use of the idea from https://github.com/zhiyhu/masonmd
#but change the input as gtf instead of variant files

#NMD is defined if satisfying all three critirea
#1) PTC is more than 50-54bp upstream of the last-exon-exon junction.
#2) targeted gene is not intronless.
#3) PTC is more than 200bp downstream of the start codon.

dictTX={} 
#txname, startcodon,cdslength,exoncount,endcodon, lastexonjunction
#final: txname, lastexondist, intronless, downstream200bp, NMD
class transcript:
	strName=strStrand="";	
	gcoord_StartCodon=None;
	gcoord_StopCodon=None;
	arrGcoordExons=[];
	arrGcoordCDS=[];
	nDist2EJ=None;
	
	def __init__(self, a_strName, a_strand):
		self.strName=a_strName;
		self.strStrand=a_strand;
		self.arrGcoordExons=[];
		self.arrGcoordCDS=[];
		self.gcoord_StartCodon=self.gcoord_StopCodon=None

	def AddCDS(self, a_start, a_end):
		self.arrGcoordCDS.append( (a_start, a_end) );

	def UpdateStartStopCodon(self):
		self.arrGcoordExons.sort();
		self.arrGcoordCDS.sort();

		nIndexLeft=nIndexRight=-1;
		arrgcoord_exon=[]
		for i in self.arrGcoordExons:
			#i is in 0-based coordinate
			#arrgcoord_exon is in 1-based coordinate
			arrgcoord_exon=arrgcoord_exon+list( range(i[0]+1, i[1]+1) )			
		
		arrgcoord_cds=[]
		for i in self.arrGcoordCDS:
			arrgcoord_cds=arrgcoord_cds+list( range(i[0]+1, i[1]+1) )
		
		nIndexLeft=arrgcoord_exon.index( arrgcoord_cds[0])
		nIndexRight=arrgcoord_exon.index( arrgcoord_cds[-1])
		gcoordStart_start=gcoordStart_end=gcoordEnd_start=gcoordEnd_end=-1;

		if self.strStrand=="+":
			gcoordStart_start=arrgcoord_exon[ nIndexLeft ]
			gcoordStart_end=arrgcoord_exon[ nIndexLeft+2 ]

			gcoordEnd_start=arrgcoord_exon[ nIndexRight+1 ]
			gcoordEnd_end=arrgcoord_exon[ nIndexRight+3 ]
		else:
			gcoordStart_start=arrgcoord_exon[ nIndexRight-2 ]
			gcoordStart_end=arrgcoord_exon[ nIndexRight ]

			gcoordEnd_start=arrgcoord_exon[ nIndexLeft-3 ]
			gcoordEnd_end=arrgcoord_exon[ nIndexLeft-1 ]

		self.UpdateStartCodon( gcoordStart_start-1, gcoordStart_end )
		self.UpdateStopCodon( gcoordEnd_start-1, gcoordEnd_end )


	def UpdateStartCodon( self, a_start, a_end ):
		self.gcoord_StartCodon=(a_start, a_end)	

	
	def UpdateStopCodon( self, a_start, a_end):
                self.gcoord_StopCodon=(a_start, a_end)
		
	def AddExon(self, a_start, a_end):
		self.arrGcoordExons.append( (a_start, a_end) );		

	def UpdateDist2EJ(self):
		if len(self.arrGcoordExons)<=1:
			return None;

		#Compared to the last EJ, what is the relative location of stop codon?
		#negative


		##	    3	      8		1-based 
		##	    2-3	      7-8	0-based

		##    EEEEEEE---------SSS     gtf	value	
		##    EEEEEEE---------E		1	>0
		##    EEEEEEE			2	=
		##    EEEE    			3	<0


		gcoordStopCodon=self.gcoord_StopCodon;
		gcoordLastEJ=None;
		nDist=None;

		if self.strStrand=="+":
			gcoordLastExon=self.arrGcoordExons[-2];
			gcoordLastEJ=(gcoordLastExon[1], gcoordLastExon[1]+1)

			nDist=gcoordStopCodon[0]-gcoordLastEJ[0]
		else:
			gcoordLastExon=self.arrGcoordExons[1];
			gcoordLastEJ=(gcoordLastExon[0], gcoordLastExon[0]+1)
		
			nDist=gcoordLastEJ[1]-gcoordStopCodon[1];
		self.nDist2EJ=nDist
			


				 
	def GetResult(self):
		nDist2SEJ=None
		nExonCount=len( self.arrGcoordExons );
		nTXLength_StopCodon=0;

		bCoding=(self.gcoord_StopCodon!=None and self.gcoord_StartCodon!=None);

		nTXLength=0;
		nCDSLength=0;
		nLeftCDS=nRigthCDS=0;
		if bCoding:
			nLeftCDS=self.gcoord_StartCodon[0] if self.strStrand=="+" else self.gcoord_StopCodon[0]+3;
			nRigthCDS=self.gcoord_StopCodon[1]-3 if self.strStrand=="+" else self.gcoord_StartCodon[1];

		for i in self.arrGcoordExons:
			nTXLength=nTXLength+i[1]-i[0];
		
			if not bCoding:
				continue;

			if i[1]<=nLeftCDS or i[0]>=nRigthCDS:
				continue;
			
			nCDSLength=nCDSLength+max(0, min(i[1], nRigthCDS)-max(i[0], nLeftCDS) );

		#Calculate nDist2SEJ
		nTXLength_StopCodon=nLastExonLen=0
		if bCoding:
			if self.strStrand=="+":
				for i in self.arrGcoordExons:
					if i[0]>=self.gcoord_StopCodon[1]:
						continue;
		
					nTXLength_StopCodon=nTXLength_StopCodon+max( min( self.gcoord_StopCodon[0], i[1] )-i[0], 0)
					nLastExonLen=self.arrGcoordExons[-1][1]-self.arrGcoordExons[-1][0]
			else:
				for i in self.arrGcoordExons:
					if i[1]<=self.gcoord_StopCodon[0]:
						continue;
					nTXLength_StopCodon=nTXLength_StopCodon+max( i[1]-max( self.gcoord_StopCodon[1], i[0] ), 0)	
					nLastExonLen=self.arrGcoordExons[0][1]-self.arrGcoordExons[0][0]
#		print("\t".join([str(nTXLength), str(nLastExonLen), str(nTXLength_StopCodon) ]))
	
		#  |-------||------------||---------||------->	 transcript (| exon junction)
		#        S                              T        stop codon
		#	  ***********************		 predicted CDS
		#  
		#  
		# Dist2SEJ : stop codon locates at upsteam of the last EJ if negative
		#	     stop codon locates at the downstream of the last EJ if positive
		
		nDist2SEJ=self.nDist2EJ if nExonCount>1 and bCoding else "NA";	
		strResult="NA"
		if bCoding:
			strResult="NotNMD"

			if nExonCount>1 and nCDSLength>200 and nDist2SEJ<-55:
				strResult="NMD"


		return "\t".join( [self.strName, "coding" if bCoding else "noncoding", str(nExonCount), str(nTXLength), str(nCDSLength), str(nDist2SEJ), strResult]);

print("Load GTF")
for gtf in tabixfile.fetch(  parser=pysam.asGTF() ):
	#loadded coordinate is 0-based 
	if gtf.feature =="gene":
		continue;

	#if gtf.transcript_id !="ENST00000648873.1":
	#	continue;

	if gtf.transcript_id not in dictTX:
		dictTX[gtf.transcript_id ]=transcript( gtf.transcript_id, gtf.strand );

	if gtf.feature == "exon":
		dictTX[gtf.transcript_id ].AddExon( gtf.start, gtf.end );
	elif gtf.feature == "CDS":
		dictTX[gtf.transcript_id ].AddCDS( gtf.start, gtf.end );

#	elif gtf.feature == "start_codon":
#		dictTX[gtf.transcript_id ].UpdateStartCodon( gtf.start, gtf.end );	
#	elif gtf.feature == "stop_codon":
#		dictTX[gtf.transcript_id ].UpdateStopCodon( gtf.start, gtf.end );
	else:
		#CDS, UTR, transcript
		continue
print("Find NMD")
for txid in dictTX:
	dictTX[ txid ].UpdateStartStopCodon();
	dictTX[ txid ].UpdateDist2EJ();


#Column	Name	Content
#1	contig	the chromosome name
#2	feature	The feature type
#3	source	The feature source
#4	start	genomic start coordinate (0-based)
#5	end	genomic end coordinate (0-based)
#6	score	feature score
#7	strand	strand
#8	frame	frame
#9	attributes	the attribute field

#

ofile=open(ofpath, "w+")
ofile.write("\t".join(["readname", "iscoding", "exoncount", "txlength", "CDSlength", "Dist2EJ", "NMDresult"])+"\n");
for tx in dictTX:
	ofile.write( dictTX[ tx ].GetResult() +"\n" );
ofile.close();

