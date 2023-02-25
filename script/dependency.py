import subprocess as sp


class LoadDependency():
	sp.Popen(['module load',"samtools/1.9"]);
	sp.Popen(['module load',"R/3.6.1"]);
	

