ipep=$1
igtf=$2
opep=$3
###ogtf=$4
min_pep=100

odir=${opep%/*}
ofname=${opep##*/}

igtfname=${igtf##*/}
ogenepred_temp=$odir/${igtfname%.gtf}".protein_in_input.genepred"

##		CDS	NMD
##transcript	O	O
##transcirpt	O	X
##transcript	X	X



#1) Extract Genepred corresponding to ipep
protein_name=$odir/${igtfname%.gtf}".protein_in_input.txt"
gtfToGenePred -genePredExt $igtf $odir/$igtfname".genepred"
awk 'NR%2==1' $ipep | awk -v FS=':' '{ gsub(/>/, "", $1); print $1; }' > $protein_name

python /script/FilterByReadName_gtf.py genepred $odir/$igtfname".genepred" $protein_name $odir/$igtfname".hasCDS_tmp.genepred" 1
python /script/annotateCDS2Genepred.py $odir/$igtfname".hasCDS_tmp.genepred" $ipep $odir/$igtfname".CDS_tmp.genepred"

#2) Output NMD prediction
oNMD=$odir/${igtfname%.gtf}".CDS.NMDprediction.txt"
genePredToGtf file $odir/$igtfname".CDS_tmp.genepred" $odir/$igtfname".CDS_tmp.gtf"

bedtools sort -i $odir/$igtfname".CDS_tmp.gtf" | bgzip -c > $odir/$igtfname".CDS_tmp.gtf.gz"
tabix -p gff $odir/$igtfname".CDS_tmp.gtf.gz"
python /script/PredictNMD.v2.py $odir/$igtfname".CDS_tmp.gtf.gz" $oNMD

#3) Extract peptide that is not NMD
awk 'NR>1 && $7=="NotNMD"{ print $1; }' $oNMD | awk -v FS='|' '{ print $1; }' > $oNMD".NMDonly.txt"
python /script/FilterByReadName_gtf.py fastq $ipep $oNMD".NMDonly.txt" ${opep%.pep}".temp.pep" 0


#4) filter by min size of aminoacid
seqtk seq -L $min_pep ${opep%.pep}".temp.pep" > $opep

#5 
ogtf=$odir/$igtfname".CDS.gtf"
python /script/annotateCDS2Genepred.py $odir/$igtfname".genepred" $opep $ogtf".genepred"
genePredToGtf file $ogtf".genepred" $ogtf

bedtools sort -i $ogtf | bgzip -c > $ogtf".gz"
tabix -p gff $ogtf.gz

rm $igtf".hasCDS_tmp.genepred.filtered" $igtf".hasCDS_tmp.genepred" $igtf".CDS_tmp.genepred" $igtf".CDS_tmp.gtf" $ogtf".genepred"
rm $oNMD".NMDonly.txt" 
rm ${opep%.pep}".temp.pep"

#



