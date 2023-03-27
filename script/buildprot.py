import os
import subprocess as sp


src_dir= os.path.dirname(os.path.realpath(__file__))

def buildprot( args ):
	print("Build proteome");
	gtf=args.gtf
	ref=args.ref
	oproteome=args.oproteome
	refproteome=args.refproteome
	thread=args.thread
	txinfo=args.txinfo
	refgtf=args.refgtf


	print(oproteome)
	process=sp.Popen(['bash', src_dir+'/prepare_proteome.sh', gtf, ref, oproteome, str(thread), txinfo, refproteome, refgtf]);
	process.wait();

	print("---Completed")









