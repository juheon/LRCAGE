IPEAK=$1
CTSSBED=$2
HALFPEAKSIZE=$3 #50 #100bp peak


#If peak size differs across technologies, precision and recall comparison is not fair
##precision: multiple number of peaks for the same TSS
##recall: closest distance from peak is sensitive to peak size

#input BED6+3

#output BED6+5
##10th column summit position
##11th colum summit read count

#Find summit
##Within peak, highest CTSS locus
ml bedtools
bedtools intersect -s -wao -a $IPEAK -b $CTSSBED | sort -k1,1 -V -k2,2n > $IPEAK".CTSSinPeak"
python2 /scratch/jmaeng/longreadCAGE/h1299/cage/script/FindSummit.v2.py $IPEAK".CTSSinPeak" $IPEAK".RepSummit"

#Extend +/- HALFPEAKSIZE to produce representative peak
awk -v OFS='\t' -v ext=$HALFPEAKSIZE -v startPos=0 '{ startPos=$10-ext; if(startPos<0){ startPos=0; } print $1,startPos,$10+ext,$4,$5,$6,$7,$8,$9,$10,$11; }' $IPEAK".RepSummit" | sort -k1,1 -V -k2,2n > ${IPEAK%.bed}".fixwd.bed"

rm $IPEAK".CTSSinPeak" $IPEAK".RepSummit"




