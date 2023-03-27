


inputdir=<input directory>
outputdir=<output directory>
scriptdir=<your LRCAGE installation directory>/script

#calling peak
docker run -v $inputdir:/inputdir -v $outputdir:/outputdir -v $scriptdir:/script lrcage \
        lrcage.py callpeak --inputlist /inputdir/<input file list> --peak /outputdir/<output prefix> --tpm 0.3 --gcap 0.35

#calling novel confident transcript
docker run -v $inputdir:/inputdir -v $outputdir:/outputdir -v $scriptdir:/script lrcage \
	lrcage.py filtertx --gtf /inputdir/<talon gtf file> --talon /inputdir/<talon tsv file> \
	--libinfo /inputdir/<library information file> --peak /outputdir/<peak file> --oprefix /outputdir/<output prefix>

#builing proteome from novel confident trancsript
docker run -v $inputdir:/inputdir -v $outputdir:/outputdir -v $scriptdir:/script lrcage \
       lrcage.py buildprot --gtf /outputdir/<novel confident trancript gtf> \
	--ref /inputdir/<reference fasta> \
	--refgtf /inputdir/<gencode gtf> \
	--txinfo /outputdir/<novel transcript information>  \
	 --refproteome /inputdir/<gencode proteome> \
	--oproteome /outputdir/<output prefix>

