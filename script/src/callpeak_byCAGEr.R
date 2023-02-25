library('CAGEr')
library('BSgenome.Hsapiens.UCSC.hg38')

##Update v3
##Use of interquantile width

##Version 5
##chr2:206159171-206162121
##changed two parameters aggregateTagClusters (maxDist to 50 from 100, excludeSignalBelowThreshold=T from F)
##Let's see whether it changes either precision/sensitivities

args=commandArgs( trailingOnly=T )
ifpath=args[1]
iIsCutoffTPM=args[2]
tagcutoff=as.numeric( args[3] )

nThread=as.numeric( args[4] );
ofpath=args[5]

options(scipen=10)

#Prepare CTSS and unG CTSS input files for peak calling
df.input=read.table(ifpath, stringsAsFactor=F, header=T, sep="\t")
df.input$samplename=apply(df.input, 1, function(a_x){ 
	strSample=paste( c(a_x[2], a_x[3], a_x[4]), collapse="_");
	print(strSample)
	return(strSample);

})

vecInputBAM=basename( df.input$input_bam )
odir=dirname( ofpath )

vecInputCTSS=stringr::str_replace( vecInputBAM, ".bam", ".bam.CTSS")
vecInputUNGCTSS=stringr::str_replace( vecInputBAM, ".bam", ".bam_unannotatedG.CTSS")

vecFName.InputCTSS=paste(odir, vecInputCTSS, sep="/")
vecFName.InputUNGCTSS=paste(odir, vecInputUNGCTSS, sep="/")


vecSampleLabels=df.input$samplename

suffix=paste0(".consclus_", tagcutoff)

cs = new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg38", inputFiles = vecFName.InputCTSS, inputFilesType = "ctss", sampleLabels = vecSampleLabels )
getCTSS(cs)
if(iIsCutoffTPM=="TRUE"){
	print("TPM peak calling")
	normalizeTagCount(cs, method="simpleTpm")
	clusterCTSS(object = cs, method="paraclu",
            threshold = tagcutoff, nrPassThreshold = 1, thresholdIsTpm = TRUE, removeSingletons = TRUE,
            keepSingletonsAbove=0.3, minStability = 2, maxLength = 100,
            reduceToNonoverlapping = TRUE, useMulticore = TRUE, nrCores = nThread )
}else{
	print("raw count peak calling")
	normalizeTagCount(cs, method="none")
	clusterCTSS(object = cs, method="paraclu",
            threshold = tagcutoff, nrPassThreshold = 1, thresholdIsTpm = FALSE, removeSingletons = TRUE,
            keepSingletonsAbove=3, minStability = 2, maxLength = 100,
            reduceToNonoverlapping = TRUE, useMulticore = TRUE, nrCores = nThread )
}

#ParaClu
cumulativeCTSSdistribution(cs, clusters="tagClusters", useMulticore = TRUE, nrCores = nThread)
quantilePositions(cs, clusters="tagClusters", qLow = 0.1, qUp = 0.9,  useMulticore = TRUE, nrCores = nThread)

tc_raw=tagClusters(cs)[[1]]
tc_inQ=tagClusters(cs, returnInterquantileWidth=T, qLow=0.1, qUp=0.9)[[1]]

tc_raw$dominant_ctss=tc_raw$dominant_ctss-1;
tc_inQ$dominant_ctss=tc_inQ$dominant_ctss-1;
tc_inQ$q_0.1=tc_inQ$q_0.1-1 #this value is 1-based coordinate so convert to 0-based

write.table( tc_raw, paste0(ofpath, ".paraclu.raw.out"), quote=F, col.names=T, row.names=F, sep="\t") ##jh
write.table( tc_inQ, paste0(ofpath, ".paraclu.interQ.out"), quote=F, col.names=T, row.names=F, sep="\t")

#Aggregate TCs within 100bp
##tag count 3 as the cutoff
##1) tagclusters with < tpmThreshold is ignored for consensus.cluster
##2) For each consensus.cluster, the peak signal is the sum of all overlapping tagclusters if excludeSignalBelowThreshold = FALSE. 
##                                                             all overlapping tagclusters with >= tpmThreshold
aggregateTagClusters(cs, tpmThreshold = tagcutoff, qLow = 0.1, qUp = 0.9, maxDist = 50,
                     excludeSignalBelowThreshold = TRUE)
                     #excludeSignalBelowThreshold = TRUE)
                     
cc_inQ=consensusClusters(cs, sample = NULL, returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
cc_inQ$clusterName=paste("cluster", cc_inQ$consensus.cluster, sep="");
cc_inQ$start=cc_inQ$start-1;
write.table( cc_inQ, paste0(ofpath, paste0(suffix, ".out") ), quote=F, col.names=T, row.names=F, sep="\t") ##jh

#BED format
##Note this score is the output of aggrateclusters
write.table( cc_inQ[, c("chr","start","end","clusterName","tpm","strand")], paste0(ofpath, paste0(suffix, ".bed") ), quote=F, col.names=F, row.names=F, sep="\t")



