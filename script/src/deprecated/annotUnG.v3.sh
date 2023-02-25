IPEAK=$1
CTSSBED=$2
UGCTSSBED=$3
gcapcutoff=$4
gcaptag=$5

ml bedtools
cut -f1-6 $IPEAK | bedtools intersect -s -wao -a stdin -b $CTSSBED | sort -k1,1 -V -k2,2n > $IPEAK".read_in_peak"
cut -f1-6 $IPEAK | bedtools intersect -s -wao -a stdin -b $UGCTSSBED | sort -k1,1 -V -k2,2n > $IPEAK".unGread_in_peak"

ml R/3.5.1 
Rscript /scratch/jmaeng/longreadCAGE/h1299/cage/script/CalcUnGRatio.v3.R $IPEAK $IPEAK".read_in_peak" $IPEAK".unGread_in_peak" ${IPEAK%.bed} $gcapcutoff $gcaptag

rm $IPEAK".read_in_peak" $IPEAK".unGread_in_peak"


