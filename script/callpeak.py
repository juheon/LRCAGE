import os
from bam2CTSS_noMAPQfilt import convertbam2CTSS
import subprocess as sp


src_dir= os.path.dirname(os.path.realpath(__file__))

def bam2CTSS(a_filelist, a_odir):
	for ifile in a_filelist:
		convertbam2CTSS( ifile, a_odir );		

	
def callpeak_cand(a_filepath, a_bCutoffTPM, a_cutoff, a_thread,  a_odir):
	process=sp.Popen(['Rscript', src_dir+"/callpeak_byCAGEr.R", a_filepath, str(a_bCutoffTPM), str(a_cutoff), str(a_thread),  a_odir ])
	process.wait();



def filter_by_gcap( a_filelist, a_odir, a_opeak,  a_gcap, a_gcap_mincount, a_half_peak_width ):
	for ifile in a_filelist:
		ifilename=ifile.rstrip("\n").split("/")[-1]
		ctss=a_odir+"/"+ifilename+".CTSS"
		ctss_unG=a_odir+"/"+ifilename+"_unannotatedG.CTSS"
		process=sp.Popen(['bash', src_dir+'/postpeakcall.v4.sh', a_opeak, ctss, ctss_unG, str(a_gcap), str( a_gcap_mincount ), str(a_half_peak_width)]);
		process.wait();	

def loadfilelist(a_idir, a_inputpath):
	arr_input=[];
	bHeader=True;
	for line in open( a_inputpath ):
		arr=line.rstrip("\n").split("\t");
		if bHeader==True:
			bHeader=False;
			continue;
		arr_input.append( a_idir +"/"+arr[0] )	
	return(arr_input);
	

def callpeak( args ):
	print("Running callpeak")
	iinputlistpath=args.inputlist
	opeak=args.peak
	thread=args.thread
	tpm=args.tpm
	readcount=args.readcount
	gcap=args.gcap
	gcap_mincount=args.gcap_mincount
	half_peak_width=args.half_peak_width

	bIsCutoffTPM=True if tpm!=None else False;	

	odir=os.path.dirname( opeak )
	sp.Popen(['mkdir', odir]);

	idir=os.path.dirname( iinputlistpath ) 
	arr_input=loadfilelist( idir, iinputlistpath )

	print("---Convert bam 2 CTSS")
	bam2CTSS( arr_input, odir );

	print("---Calling peaks")
	callpeak_cand( iinputlistpath, bIsCutoffTPM, tpm if bIsCutoffTPM else readcount, thread, opeak);

	print("---Filtering by G-cap")
	opeak_raw=opeak+".consclus_"+str(tpm)+".bed"					
	filter_by_gcap( arr_input, odir, opeak_raw, gcap, gcap_mincount, half_peak_width )
#
	print("---Completed")
	



