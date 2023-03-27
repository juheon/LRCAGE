gtf=$1
talon=$2
libinfo=$3
mincount=$4
minpeakratio=$6
peak=$5
outprefix=$7


odir=${outprefix%/*}
gtfname=${gtf##*/}

#1. Calculate the number of read per transcript per dataset. Annotate whether the number of read passes the minimum readcount
Rscript /script/Talon2TXinfo.R $talon $mincount $outprefix.readcount_pertx.txt

#2. Filtering by peak
tolerance=250
## 1) Count the number of read connecting peak and transcript
awk -v OFS='\t' -v tolerance=$tolerance -v start=-1 -v end=-1 'NR>1{ 
	start=$5-tolerance; end=$5+tolerance;
	if(start<0){ start=0; } 
	print $4,start,end,$1,0,$7,$2,$13;}' $talon | sort -k1,1 -V -k2,2n > $outprefix.read_TSS.bed
bedtools intersect -wo -s -a <(cut -f1-6 $peak | sort -k1,1 -V -k2,2n) -b $outprefix.read_TSS.bed > $outprefix.peak_with_read_TSS.bed
Rscript /script/count_read_per_peaktx.R $outprefix.peak_with_read_TSS.bed $outprefix.readcount_per_txpeak.txt
rm $outprefix.read_TSS.bed $outprefix.peak_with_read_TSS.bed



## 2) Find the best peak for each transcript
## Calculate TPM
Rscript /script/find_best_peak.R $outprefix.readcount_pertx.txt $outprefix.readcount_per_txpeak.txt $libinfo $minpeakratio $outprefix.txinfo.txt
rm $outprefix.readcount_pertx.txt $outprefix.readcount_per_txpeak.txt


#3. Extract confident gtf
#gtfToGenePred -genePredExt $gtf $outprefix.genepred
python /script/update_genepred.py $outprefix.genepred $outprefix.txinfo.txt $outprefix True
genePredToGtf file -source="conftx" $outprefix.conf.genepred $outprefix.conf.gtf
bedtools sort -i $outprefix.conf.gtf | bgzip -c > $outprefix.conf.gtf.gz; tabix -p gff $outprefix.conf.gtf.gz
rm $outprefix.conf.genepred $outprefix.conf.gtf $outprefix.genepred









 



 





