


inputdir=<input directory>
outputdir=<output directory>
scriptdir=<LRCAGE installation directory>/script

##Calling peak
docker run -v $inputdir:/inputdir -v $outputdir:/outputdir -v $scriptdir:/script lrcage \
        lrcage.py callpeak --inputlist /inputdir/inputlist.txt --peak /outputdir/<output prefix> --tpm 0.3 --gcap 0.35

##Profiling novel confident transcripts from LRCAGE data
docker run -v $inputdir:/inputdir -v $outputdir:/outputdir -v $scriptdir:/script lrcage \
	lrcage.py filtertx --gtf /outputdir/<TALON gtf> --talon /outputdir/<TALON tsv> \
	--libinfo /inputdir/<library information txt> --peak /outputdir/<peak bed> --oprefix /outputdir/<output prefix>

##Profiling proteome from novel confident transcripts
docker run -v $inputdir:/inputdir -v $outputdir:/outputdir -v $scriptdir:/script lrcage \
       lrcage.py buildprot --gtf /outputdir/<novel confident transcript gtf> \
	--ref /inputdir/<reference genome fasta> \
	--refgtf /inputdir/<gencode gtf> \
	--txinfo /outputdir/<transcript information txt> \
	 --refproteome /inputdir/<gencode proteome fasta> \
	--oproteome /outputdir/<output prefix>

