# programmer : Nakul, adapted from code written by Daofeng
# Ju Heon modified to use more lenient criteria of G cap filter. The original version only use "G", "GG", "GGG" as the G cap filtering but I changed this to "G", "G*", "G**".
# usage: This is a quick method to generate a CTSS file and a unannotatedG CTSS file
# Other update 10212020
## MAPQ filter is already applied in input bam file

import subprocess as sp
import re

BAM_FPAIRED = 1 
BAM_FPROPER_PAIR = 2 
BAM_FUNMAP = 4
BAM_FMUNMAP = 8 
BAM_FREVERSE = 16
BAM_FMREVERSE = 32 
BAM_FREAD1 = 64
BAM_FREAD2 = 128 
BAM_FSECONDARY = 256 

#HWI-ST841:229:C2EPTACXX:7:1102:2740:73390       163     chr5    138386709       255     27M69988N74M    =       138463484       147287  CGTGGGATCCTGCAGGGACAGCCTGCCCTGGCTGAAGGGCCTGCCACTCATGCGTCGGGTGGAACACCTCCAGGA
#HWI-ST841:229:C2EPTACXX:7:1102:2740:73390       83      chr5    138463484       255     59M70415N38M4S  =       138386709       -147287 GCATGCCCAGAGGAGCCATCCTAGATGAAGGCAGGCTCTGGGGAGCCATAGTCAGGGACCTTCGGCCTCCGGCAG
def bam_calend(start, cigar):
    """
    >>> a = '1S22M98N25M'
    >>> a
    '1S22M98N25M'
    >>> import re
    >>> re.split('\d+', a)
    ['', 'S', 'M', 'N', 'M']
    >>> re.split('\D+', a)
    ['1', '22', '98', '25', '']
    """
    ops = re.split('\d+', cigar)
    cnts = re.split('\D+', cigar)
    ops = ops[1:]
    cnts = cnts[:-1]
    toadd = 0
    for i in range(len(ops)):
        if ops[i] != 'I' and ops[i] != 'S':
            toadd += int(cnts[i])
    return start + toadd -1 #jmaeng: For single-end seq, TLEN field is always zero and bam_calend should be used to calcualte the template size. toadd calculcates using the CIGAR string. Since the coordinate is 1-based coordinate, the end position should be start+add-1 not start+add


def convertbam2CTSS( a_inf, a_odir):
    ud = {} #for unique, key: position, chr1:5000, value: tag count
    udG = {} #for unique reads from unannotated Gs, key: position, chr1:5000, value: tag count
    udGName = [] #for the names of the unannotated G reads so they can be extracted from the file

    uc = 0
    c = 0
    inf = a_inf
    odir = a_odir #jmaeng
    infk = inf.split('/')[-1] #get the filename but do not put everything in the foreign folder
    p = sp.Popen(['samtools','view', inf], stdout=sp.PIPE)
    prid = ''
    for line in p.stdout:
        strLine=line.decode('ASCII');

        if strLine.startswith('@'): continue
        #print line.strip()
        #continue
        t = strLine.strip().split('\t')
        rid = t[0]
        flag = int(t[1])
        chrom = t[2]
        cigar = t[5]
        mapQ = int(t[4])
        isize = int(t[8])
        seqQ = t[9]
        
        #if mapQ == 255: #Only unique hits are used. #Jh: This is STAR-specific definition which may not hold for other aligners. By SAM format document, 255 is forbidden value. 
	#Minimap uses 255 for the missing alignment. jmaeng
	#if mapQ != 255:
        if flag & BAM_FREVERSE:
            strand = '-'
        else:
            strand = '+'
            
        if flag & BAM_FPAIRED:
            if flag & BAM_FREAD1:
                if rid != prid:
                    c += 1
                if strand == '+':
                    start = int(t[3])
                else:
                    start = int(t[7])
            else:
                continue
        else:
            if rid != prid:
                c += 1
            start = int(t[3])
        if isize != 0:
            if strand == '+':
                end = start + isize -1
            else:
                end = start - isize -1
        else:
            end = bam_calend(start, cigar)
        if strand == '+':
            pos = '{}:{}:{}'.format(chrom, start, strand)
        else:
            pos = '{}:{}:{}'.format(chrom, end, strand)
        if rid != prid:
            uc += 1
            if pos not in ud:
                ud[pos] = 1
            else:
                ud[pos] += 1
            
            #This is to see if there are unannotated G's in forward or reverse (C's) direction
            if strand == '+':
                if cigar[0:2] == '1S' and seqQ[0] == 'G':
                    if pos not in udG:
                        udG[pos] = 1
                    else:
                        udG[pos] += 1
                if cigar[0:2] == '2S' and seqQ[0] == 'G': #seqQ[0:2] == 'GG':
                    if pos not in udG:
                        udG[pos] = 1
                    else:
                        udG[pos] += 1
                if cigar[0:2] == '3S' and seqQ[0] == 'G': #seqQ[0:3] == 'GGG':
                    if pos not in udG:
                        udG[pos] = 1
                    else:
                        udG[pos] += 1
            else:
                if cigar[-2:] == '1S' and seqQ[-1] == 'C':
                    if pos not in udG:
                        udG[pos] = 1
                    else:
                        udG[pos] += 1
                if cigar[-2:] == '2S' and seqQ[-1] == 'C': #seqQ[-2:] == 'CC':
                    if pos not in udG:
                        udG[pos] = 1
                    else:
                        udG[pos] += 1
                if cigar[-2:] == '3S' and seqQ[-1] == 'C': #seqQ[-3:] == 'CCC':
                    if pos not in udG:
                        udG[pos] = 1
                    else:
                        udG[pos] += 1
                    
        ##else:
        ##    continue
        if c % 1000000 == 0:
            print(" ".join(['processed', str(c), 'tags']));
        prid = rid


    with open(odir+"/"+'{}.CTSS'.format(infk), 'w') as fbedu:
        for i in ud:
            t = i.split(':')
            if t[2] == '+':
                fbedu.write('{}\t{}\t+\t{}\n'.format(t[0], t[1],  ud[i]))
            else:
                fbedu.write('{}\t{}\t-\t{}\n'.format(t[0], t[1], ud[i]))
                
    with open(odir+"/"+'{}_unannotatedG.CTSS'.format(infk), 'w') as fbedu:

        for i in udG:
            t = i.split(':')
            if t[2] == '+':
                fbedu.write('{}\t{}\t+\t{}\n'.format(t[0], t[1],  udG[i]))
            else:
                fbedu.write('{}\t{}\t-\t{}\n'.format(t[0], t[1], udG[i]))


if __name__=="__main__":
    main()


