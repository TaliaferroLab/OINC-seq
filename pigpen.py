#Pipeline for Identification of Guanosine Positions Erroneously Notated
#PIGPEN

import argparse
import subprocess
import os
from snps import getSNPs, recordSNPs
from filterbam import intersectreads, filterbam
from getmismatches import iteratereads_pairedend, getmismatches
from assignreads import getReadOverlaps, processOverlaps
from conversionsPerGene import getPerGene, writeConvsPerGene

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'pigpen for quantifying OINC-seq data')
    parser.add_argument('--bam', type = str, help = 'Aligned reads (ideally STAR uniquely aligned reads) to quantify',
                        required = True)
    parser.add_argument('--controlBams', type = str, help = 'Comma separated list of alignments from control samples (i.e. those where no *induced* conversions are expected. Required if SNPs are to be considered.')
    parser.add_argument('--genomeFasta', type = str, help = 'Genome sequence in fasta format. Required if SNPs are to be considered.')
    parser.add_argument('--geneBed', type = str, help = 'Bed file of genomic regions to quantify. Fourth field must be gene ID.')
    parser.add_argument('--chromsizes', type = str, help = 'Tab-delimited file of chromosomes in the order the appear in the bed/bam and their sizes. Can be made with cut -f 1,2 genome.fa.fai')
    parser.add_argument('--output', type = str, help = 'Output file of conversion rates for each gene.')
    parser.add_argument('--nproc', type = int, help = 'Number of processors to use. Default is 1.', default = 1)
    parser.add_argument('--useSNPs', action = 'store_true', help = 'Consider SNPs?')
    parser.add_argument('--SNPcoverage', type = int, help = 'Minimum coverage to call SNPs. Default = 20', default = 20)
    parser.add_argument('--SNPfreq', type = float, help = 'Minimum variant frequency to call SNPs. Default = 0.02', default = 0.02)
    parser.add_argument('--onlyConsiderOverlap', action = 'store_true', help = 'Only consider conversions seen in both reads of a read pair?')
    parser.add_argument('--requireMultipleConv', action = 'store_true', help = 'Only consider conversions seen in reads with multiple G->C + G->T conversions?')
    args = parser.parse_args()

    #Make index for bam if there isn't one already
    bamindex = args.bam + '.bai'
    if not os.path.exists(bamindex):
        indexCMD = 'samtools index ' + args.bam
        index = subprocess.Popen(indexCMD, shell = True)
        index.wait()
    
    #Make vcf file for snps
    if args.useSNPs:
        controlbams = args.controlBams.split(',')
        
        #Make index for each control bam if there isn't one already
        for bam in controlbams:
            bamindex = bam + '.bai'
            if not os.path.exists(bamindex):
                indexCMD = 'samtools index ' + bam
                index = subprocess.Popen(indexCMD, shell = True)
                index.wait()

        vcfFileNames = getSNPs(controlbams, args.genomeFasta, args.SNPcoverage, args.SNPfreq)
        snps = recordSNPs('merged.vcf')
    
    elif not args.useSNPs:
        snps = None

    #Filter bam for reads contained within entries in geneBed
    #This will reduce the amount of time it takes to find conversions
    print('Filtering bam for reads contained within regions of interest...')
    intersectreads(args.bam, args.geneBed, args.chromsizes)
    filteredbam = filterbam(args.bam)

    #Identify conversions
    if args.nproc == 1:
        convs, readcounter = iteratereads_pairedend(filteredbam, args.onlyConsiderOverlap, snps, args.requireMultipleConv, 'high')
    elif args.nproc > 1:
        convs = getmismatches(filteredbam, args.onlyConsiderOverlap, snps, args.requireMultipleConv, args.nproc)


    #Assign reads to genes
    print('Assigning reads to genes...')
    overlaps, numpairs = getReadOverlaps(filteredbam, args.geneBed, args.chromsizes)
    read2gene = processOverlaps(overlaps, numpairs)

    #Calculate number of conversions per gene
    numreadspergene, convsPerGene = getPerGene(convs, read2gene)
    writeConvsPerGene(numreadspergene, convsPerGene, args.output)





        

