gtf=$1
ref=$2
oprefix=$3
thread=$4
txinfo=$5
refproteomefasta=$6
refgtf=$7


odir=${oprefix%/*}
mkdir $odir




ml kentUCSC bedtools samtools R/3.6.1 python3
#1. Convert gtf to fasta
gtfname=${gtf##*/}
#gtfToGenePred -genePredExt $gtf $odir/$gtfname.genepred
#genePredToBed $odir/$gtfname.genepred $odir/$gtfname.bed12
#bedtools getfasta -name -s -split -fi $ref -bed $odir/$gtfname.bed12 > $odir/$gtfname.fasta

#2. Removal of fasta with non-ATCG
#python /bar/jmaeng/pfiles/LRCAGE/script/src/FindNonATCG.py $odir/$gtfname.fasta $odir/$gtfname.removeN.fasta
#rm $odir/$gtfname.genepred $odir/$gtfname.bed12 $odir/$gtfname.fasta

#3. Run ORF
#dumb_predict.py --cpus $thread $odir/$gtfname.removeN.fasta $odir/$gtfname.removeN.longORF_cand.fasta
#python /bar/jmaeng/pfiles/LRCAGE/script/src/FilterORF_angel.py $odir/$gtfname.removeN.longORF_cand.fasta.final.pep 1 1 $odir/$gtfname.removeN.longORF.final.pep


#4. Discard NMD
#bash /bar/jmaeng/pfiles/LRCAGE/script/src/predictNMD.v4.sh $odir/$gtfname.removeN.longORF.final.pep $gtf $odir/$gtfname.removeN.longORF.final.noNMD.pep

#5. Deduplicate proteins
iprotein=$odir/$gtfname.removeN.longORF.final.noNMD.pep
#python /bar/jmaeng/pfiles/LRCAGE/script/src/prepare_proteome.v2.py $iprotein $txinfo $oprefix
#samtools faidx $oprefix".proteome.fasta"
#rm $gtf.genepred  $gtf.bed12 ${oprefix%.fasta}.removeN.fasta ${oprefix%.fasta}.removeN.longORF_cand.fasta

#6. Compare protein against GENCODE proteome
iproteomefasta=$oprefix".proteome.fasta"
iproteometable=$oprefix".proteome_table.txt"
#python /bar/jmaeng/pfiles/LRCAGE/script/src/compare_frame_score.v8.py $iproteomefasta $iproteometable $refproteomefasta $oprefix

gtf_cds=$odir/$gtfname".CDS.gtf.gz"
python /bar/jmaeng/pfiles/LRCAGE/script/src/filter_spurious.v4.py $oprefix".framecall_align.txt" $gtf_cds $refgtf $oprefix".framecall_align.filt.txt"



#6. Extract protein coordinate of unannotated peptide
python /bar/jmaeng/pfiles/LRCAGE/script/src/find_eventpcoord_againstRef.v3.py $oprefix".framecall_align.filt.txt" $oprefix".framecall_align.filt"



#d 7. Combine TALON output with protein information
txinfo_CDS=${txinfo%.txt}".CDS.txt"
Rscript /bar/jmaeng/pfiles/LRCAGE/script/src/combine_protein.R $txinfo $iproteometable $oprefix".framecall_align.filt.stat.txt" $txinfo_CDS





