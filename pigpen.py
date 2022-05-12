#Pipeline for Identification of Guanosine Positions Erroneously Notated
#PIGPEN

import argparse
import subprocess
import os
import sys
from snps import getSNPs, recordSNPs
from filterbam import intersectreads, filterbam, intersectreads_multiprocess
from getmismatches import iteratereads_pairedend, getmismatches
from assignreads import getReadOverlaps, processOverlaps
from conversionsPerGene import getPerGene, writeConvsPerGene
import pickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='                    ,-,-----,\n    PIGPEN     **** \\ \\ ),)`-\'\n              <`--\'> \\ \\` \n              /. . `-----,\n    OINC! >  (\'\')  ,      @~\n              `-._,  ___  /\n-|-|-|-|-|-|-|-| (( /  (( / -|-|-| \n|-|-|-|-|-|-|-|- \'\'\'   \'\'\' -|-|-|-\n-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|\n\n   Pipeline for Identification \n      Of Guanosine Positions\n       Erroneously Notated', formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bam', type = str, help = 'Aligned reads (ideally STAR uniquely aligned reads) to quantify', required = True)
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
    parser.add_argument('--use_g_t', action = 'store_true', help = 'Consider G->T conversions?')
    parser.add_argument('--use_g_c', action = 'store_true', help = 'Consider G->C conversions?')
    parser.add_argument('--nConv', type = int, help = 'Minimum number of required G->T and/or G->C conversions in a read pair in order for conversions to be counted. Default is 1.', default = 1)
    args = parser.parse_args()

    #We have to be either looking for G->T or G->C, if not both
    if not args.use_g_t and not args.use_g_c:
        print('We have to either be looking for G->T or G->C, if not both! Add argument --use_g_t and/or --use_g_c.')
        sys.exit()

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
    if args.nproc == 1:
        intersectreads(args.bam, args.geneBed, args.chromsizes)
        filteredbam = filterbam(args.bam, args.nproc)
    elif args.nproc > 1:
        intersectreads_multiprocess(args.bam, args.geneBed, args.chromsizes, args.nproc)
        filteredbam = filterbam(args.bam, args.nproc)

    #Identify conversions
    if args.nproc == 1:
        convs, readcounter = iteratereads_pairedend(filteredbam, args.onlyConsiderOverlap, args.use_g_t, args.use_g_c, snps, args.nConv, 'high')
    elif args.nproc > 1:
        convs = getmismatches(filteredbam, args.onlyConsiderOverlap, snps, args.nConv, args.nproc, args.use_g_t, args.use_g_c)


    #Assign reads to genes
    print('Assigning reads to genes...')
    overlaps, numpairs = getReadOverlaps(filteredbam, args.geneBed, args.chromsizes)
    read2gene = processOverlaps(overlaps, numpairs)
    os.remove(filteredbam)
    os.remove(filteredbam + '.bai')

    #Calculate number of conversions per gene
    numreadspergene, convsPerGene = getPerGene(convs, read2gene)
    writeConvsPerGene(numreadspergene, convsPerGene, args.output, args.use_g_t, args.use_g_c)





        

