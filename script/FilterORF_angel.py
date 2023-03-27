import sys
import pysam

ifasta=sys.argv[1];
ofasta=sys.argv[4];
bRemoveStop=int(sys.argv[2]);	#1: remove
bFilterComplete=int(sys.argv[3])#1: filter

ofile=open( ofasta, "w+");
for read in pysam.FastxFile( ifasta ):
	arr=read.comment.split(" ")	
	strType=arr[0].lstrip("type:");

	if bFilterComplete==1 and strType!="dumb-complete":
		continue;

	if bRemoveStop==1:
		read.sequence=read.sequence.rstrip("*");
	ofile.write( str( read )+"\n");

ofile.close()




