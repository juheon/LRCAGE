import os
import subprocess as sp


src_dir= os.path.dirname(os.path.realpath(__file__))

def filtertx( args ):
	print("Calling confident transcripts")
	gtf=args.gtf
	talon=args.talon
	libinfo=args.libinfo;
	peak=args.peak
	mincount=args.mincount
	minpeakratio=args.peakratio
	outprefix=args.oprefix
	
	if peak is None:
		peak="None";

	process=sp.Popen(['bash', src_dir+'/filter_transcript.sh', gtf, talon, libinfo, str(mincount), peak, str(minpeakratio), outprefix]);
	process.wait();
	print("---Completed")









