PARACLUBED=$1
CTSSBED=$2
CTSS_UG_BED=$3
fGCAP=$4
nGCAPtag=$5
HALFPEAKWIDTH=$6

##Filter by peak intensity and by gcap ratio
cut -f1-6 $PARACLUBED | bedtools intersect -s -wao -a stdin -b <(awk -v OFS='\t' '{ print $1,$2-1,$2,"CTSS",$4,$3;}' $CTSSBED | sort -k1,1 -V -k2,2n) > $PARACLUBED".read_in_peak"
cut -f1-6 $PARACLUBED | bedtools intersect -s -wao -a stdin -b <(awk -v OFS='\t' '{ print $1,$2-1,$2,"CTSS",$4,$3;}' $CTSS_UG_BED | sort -k1,1 -V -k2,2n) > $PARACLUBED".unGread_in_peak"
Rscript /script/CalcUnGRatio.v3.R $PARACLUBED $PARACLUBED".read_in_peak" $PARACLUBED".unGread_in_peak" ${PARACLUBED%.bed} $fGCAP $nGCAPtag
rm $PARACLUBED".read_in_peak" $PARACLUBED".unGread_in_peak"

##Make fixed-width: Final peak will be 2*HALFPEAKWIDTH long
#input BED6+3

#output BED6+5
##10th column summit position
##11th colum summit read count

#Find summit
##Within peak, highest CTSS locus
bedtools intersect -s -wao -a ${PARACLUBED%.bed}".unG_"$fGCAP".bed" -b <(awk -v OFS='\t' '{ print $1,$2-1,$2,"CTSS",$4,$3;}' $CTSSBED | sort -k1,1 -V -k2,2n) | sort -k1,1 -V -k2,2n > ${PARACLUBED%.bed}".unG_"$fGCAP".CTSSinPeak"
python /script/FindSummit.v2.py ${PARACLUBED%.bed}".unG_"$fGCAP".CTSSinPeak" ${PARACLUBED%.bed}".unG_"$fGCAP".RepSummit"


#Extend +/- HALFPEAKSIZE to produce representative peak
cat <(echo -e '#chrom\tstart\tend\tconsensus.cluster\treadcnt\tstrand\treadcnt\treadcnt_unG\tunGratio\tsummit\treadcnt_summit') <(awk -v OFS='\t' -v ext=$HALFPEAKSIZE -v startPos=0 '{ startPos=$10-ext; if(startPos<0){ startPos=0; } print $1,startPos,$10+ext,$4,$5,$6,$7,$8,$9,$10,$11; }' ${PARACLUBED%.bed}".unG_"$fGCAP".RepSummit" | sort -k1,1 -V -k2,2n)  > ${PARACLUBED%.bed}".unG_"$fGCAP".fixwd.bed"
rm ${PARACLUBED%.bed}".unG_"$fGCAP".RepSummit" ${PARACLUBED%.bed}".unG_"$fGCAP".CTSSinPeak"
rm ${PARACLUBED%.bed}".unG_"$fGCAP".bed" ${PARACLUBED%.bed}".unG.bed"






