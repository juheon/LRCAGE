import sys
import pysam

ifasta=sys.argv[1];
ofasta=sys.argv[2];

ofile=open( ofasta, "w+");
for read in pysam.FastxFile( ifasta ):
	if read.sequence.count("N")>0 or read.sequence.count("n")>0:
		continue;

	ofile.write( str(read)+"\n");

ofile.close();





