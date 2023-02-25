

idir=/scratch/jmaeng/longreadCAGE/h1299/LRCAGE/test
#../script/lrcage.py callpeak --inputlist $idir/inputlist.bamfiles.simpleinput.txt --peak $idir/peak/nanoCAGE --tpm 0.3 --gcap 0.35

talon=`pwd`/LRCAGE_WT_DACSB_gtf_talon_read_annot.tsv
libinfo=`pwd`/filtertx/libinfo.txt
gtf=`pwd`/LRCAGE_WT_DACSB_gtf.TALON_talon.gtf
peak=/scratch/jmaeng/longreadCAGE/h1299/pacbio/downsample.v2/consensus_peak/WT_DACSB_consenssuspeak.consclus_0.9.unG_0.35.bed
#../script/lrcage.py filtertx --gtf $gtf --talon $talon --libinfo $libinfo --peak $peak --oprefix `pwd`/filtertx/LRCAGE.DMSO_DACSB


gtf=$idir/filtertx/LRCAGE.DMSO_DACSB.conf_nc.gtf.gz
ref=~/genomes/hg38_noalt_maskPAR/hg38_noalt_maskPAR.fa
refgtf=/scratch/jmaeng/longreadCAGE/h1299/annotation/gencode.v29.basic.annotation.gtf.gz
txinfo=$idir/filtertx/LRCAGE.DMSO_DACSB.txinfo.txt
outproteome=$idir/proteome/LRCAGE.DMSO_DACSB
refproteome=/scratch/jmaeng/longreadCAGE/h1299/mass_spec/custom_proteome/GENCODE_basic_proteome.wGeneID.fasta
./lrcage.py buildprot --gtf $gtf --ref $ref --txinfo $txinfo --oproteome $outproteome --refproteome $refproteome --refgtf $refgtf


