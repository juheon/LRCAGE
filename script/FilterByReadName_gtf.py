import sys
import pysam

inputtype=sys.argv[1];
ifasta=sys.argv[2];
iconf=sys.argv[3];
ofasta=sys.argv[4];
bOutFiltered=int(sys.argv[5]);	#1

if bOutFiltered==1:
	ofile2=open( ofasta+".filtered", "w+" );

dictConf={}
for line in open( iconf ):
	strReadName=line.rstrip("\n").split("\t")[0]
	dictConf[ strReadName ]={strReadName}

if inputtype=="fastq":
	ofile=open( ofasta, "w+");
	for read in pysam.FastxFile( ifasta ):
		arr=read.name.split(":")
		strReadName=arr[0]
	
		if strReadName.count("_")>1:
			#For 3frame
			strReadName=strReadName.split("_")[0]
	
		if strReadName not in dictConf:
			if bOutFiltered==1:
				ofile2.write( str( read )+"\n");
		else:
			ofile.write( str( read )+"\n");
	ofile.close()
elif inputtype=="gtf":
	ofile=open( ofasta, "w+");
	gtffile=pysam.TabixFile( ifasta, parser=pysam.asGTF() );
	for entry in gtffile.fetch():
		if entry.feature == "gene":
			continue

		if entry.transcript_id not in dictConf:
			if bOutFiltered==1:
				ofile2.write( str( entry )+"\n");		
		else:
			ofile.write( str( entry )+"\n" );
	ofile.close()
elif inputtype=="genepred":
	ofile=open( ofasta, "w+");
	for line in open( ifasta ):
		arr=line.rstrip("\n").split("\t");
		if arr[0] not in dictConf:
			if bOutFiltered==1:
				ofile2.write( "\t".join( arr )+"\n"); 
		else:
			ofile.write( "\t".join( arr )+"\n");
	ofile.close();
else:
	print("invalid input type");


if bOutFiltered==1:
	ofile2.close();



