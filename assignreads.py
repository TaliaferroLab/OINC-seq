#Take a bam (ideally uniquely aligned and pre-filtered for alignments that fall within 
#regions of interest) and for each read, assign it to a region in a bam file.
#If the read overlaps with >1 gene, assign it to the gene it overlaps the most with.

import pybedtools
import os
import sys
import pickle
import subprocess

def getReadOverlaps(bam, bed, chrsort):
    #Given a bam and a bed of regions to overlap with,
    #for each read, print the regions in the bed it overlaps with
    #as well as the amount of overlap.
    #This is likely going to be slow.
    #Both the bam and the bed must be coordinate sorted

    #chrsort is a tab-delimited file of chromosomes in the order the appear in the bed/bam and their sizes
    #can be made with cut -f 1,2 genome.fa.fai

    #Get number of read pairs in this file (assuming this file only contains uniquely mapped reads)
    numreads = subprocess.check_output(['samtools', 'view', '-c', bam])
    numpairs = int(int(numreads) / 2) #paired end reads

    bam = pybedtools.BedTool(bam)
    bed = pybedtools.BedTool(bed)

    print('Intersecting {0} read pairs...'.format(numpairs))

    intersections = bam.intersect(bed, sorted = True, g = chrsort, bed = True, split = True, wo = True)

    overlaps = {} #{readid : {txid : overlaplength}}
    for i in intersections:
        i = str(i).strip().split('\t')
        readid = i[3].split('/')[0]
        #bedtools reports intersections for each read of a mate pair separately
        #we are going to combine their intersection lengths
        txid = i[15]
        overlaplength = int(i[18])
        if readid not in overlaps:
            overlaps[readid] = {}
        
        if txid not in overlaps[readid]:
            overlaps[readid][txid] = overlaplength
        elif txid in overlaps[readid]:
            overlaps[readid][txid] += overlaplength
    
    return overlaps, numpairs

def processOverlaps(overlaps, numpairs):
    #For each read, get the transcript with which it has the most overlap
    read2gene = {} #{readid : ensembl_gene_id}
    for read in overlaps:
        txs = overlaps[read]
        maxtx = max(txs, key = txs.get)
        overlaplength = txs[maxtx] #can implement minimum overlap here
        gene = maxtx.split('_')[0]
        read2gene[read] = gene

    frac_readpairs_with_gene = round((len(read2gene) / numpairs) * 100, 2)
    print('Found genes for {0} readpairs ({1}%).'.format(len(read2gene), frac_readpairs_with_gene))

    return read2gene



if __name__ == '__main__':
    overlaps, numpairs = getReadOverlaps(sys.argv[1], sys.argv[2], sys.argv[3])
    read2gene = processOverlaps(overlaps, numpairs)

    with open('read2gene.pkl', 'wb') as outfh:
        pickle.dump(read2gene, outfh)


