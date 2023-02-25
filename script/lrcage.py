#!/opt/apps/python3/bin/python3

import argparse


#Load src files
import os 
src_dir= os.path.dirname(os.path.realpath(__file__))+"/src"

import sys
sys.path.insert(0, src_dir);
from callpeak import *
from buildprot import * 
from filtertx import *
from annotpep import *

def range_float(a_min, a_max):

	def is_float_valid( a_x ):
		try:
			fInput=float(a_x);	
		except ValueError:
			raise argparse.ArgumentTypeError("Only float value is allowed")	

		if fInput < a_min or fInput > a_max:
			raise argparse.ArgumentTypeError("Out of range");

		return fInput;

	return is_float_valid



if __name__ == "__main__":
	parser=argparse.ArgumentParser(prog="LRCAGE", description="Analysis pipeline for long-read CAGE data", epilog="");
	parser.add_argument('-ver', action='version', version='%(prog)s 1.0');
	subparsers = parser.add_subparsers(help='Description')
	
	parser_1 = subparsers.add_parser('callpeak', help='to call peak')
	parser_4 =subparsers.add_parser('filtertx', help='to retain a list of confident transcripts');
	parser_3 =subparsers.add_parser('buildprot', help='to create a proteome database');
	parser_2 = subparsers.add_parser('annotpep', help='to annotate peptide')	


	#peak calling
	parser_1.add_argument('--inputlist', help='list of input bam files', required=True)	
	parser_1.add_argument('--peak', help='output peak file name', required=True)

	group_peakcutoff=parser_1.add_mutually_exclusive_group(required=True)
	group_peakcutoff.add_argument('--tpm', type=float, help='minimum TPM per peak')
	group_peakcutoff.add_argument('--readcount', type=int, help='minimum read count per peak')

	parser_1.add_argument('--gcap', type=range_float(0, 1), help='minimum G-cap ratio', default=0.1);
	parser_1.add_argument('--gcap_mincount', type=int, help='minimum number of soft-clipped G reads', default=2)
	parser_1.add_argument('--half_peak_width', type=int, help='half peak size', default=50)

	parser_1.add_argument('--thread', type=int, default=4, help='number of threads' );
	parser_1.set_defaults(func=callpeak)
	

	#filter peptide
	parser_2.add_argument('--peptide', help='input peptide.txt file from MaxQuant')
	parser_2.add_argument('--ref', help='reference proteome to define canonical peptides')
	parser_2.add_argument('--protinfo', help="protain class information");
	parser_2.add_argument('--opeptide', help="output file name" )
	parser_2.set_defaults(func=annotpep)


	#proteome asssembly
	parser_3.add_argument('--gtf', help="input gtf file", required=True)
	parser_3.add_argument('--ref', help="reference genome fasta", required=True)
	parser_3.add_argument('--txinfo', help="transcript information", required=False)
	parser_3.add_argument('--thread', type=int, default=4, help='number of threads' );
	parser_3.add_argument('--oproteome', help="output proteome", required=True);
	parser_3.add_argument('--refproteome', help="reference proteome", required=True);
	parser_3.add_argument('--refgtf', help="reference gtf", required=True);
	parser_3.set_defaults(func=buildprot)	

	
	#filter transcript
	parser_4.add_argument('--gtf', help="input gtf file", required=True);
	parser_4.add_argument('--talon', help="input TALON.tsv file", required=True);
	parser_4.add_argument('--libinfo', help="library size information");
	parser_4.add_argument('--mincount', help="minimum count to define confident transcripts", default=3)
	parser_4.add_argument('--peak', help="peaks used to retain transcripts with complete 5' ends", required=False);
	parser_4.add_argument('--peakratio', help="minimum fraction of reads for peak-transcript pair per trancsript", required=False, default=0.3);
	parser_4.add_argument('--oprefix', help="prefix for output files", required=True);
	parser_4.set_defaults(func=filtertx)



	try:
		args = parser.parse_args()
		if hasattr( args, 'func' ):
			args.func( args );
		else:
			parser.parse_args(['--help'])

	except argparse.ArgumentError:
		print("ERR")	



