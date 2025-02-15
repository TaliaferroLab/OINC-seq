#Pipeline for Identification of Guanosine Positions Erroneously Notated
#PIGPEN

import argparse
import subprocess
import os
import sys
from snps import getSNPs, recordSNPs
from maskpositions import readmaskbed
from getmismatches import iteratereads_pairedend, getmismatches
from assignreads_salmon_ensembl import getpostmasterassignments, assigntotxs, collapsetogene, readspergene, writeOutput
from assignreads import getReadOverlaps, processOverlaps
from conversionsPerGene import getPerGene, writeConvsPerGene

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='                    ,-,-----,\n    PIGPEN     **** \\ \\ ),)`-\'\n              <`--\'> \\ \\` \n              /. . `-----,\n    OINC! >  (\'\')  ,      @~\n              `-._,  ___  /\n-|-|-|-|-|-|-|-| (( /  (( / -|-|-| \n|-|-|-|-|-|-|-|- \'\'\'   \'\'\' -|-|-|-\n-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|\n\n   Pipeline for Identification \n      Of Guanosine Positions\n       Erroneously Notated', formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--datatype', type = str, choices = ['single', 'paired'], required = True, help = 'Single end or paired end data?')
    parser.add_argument('--samplenames', type = str, help = 'Comma separated list of samples to quantify.', required = True)
    parser.add_argument('--controlsamples', type = str, help = 'Comma separated list of control samples (i.e. those where no *induced* conversions are expected). May be a subset of samplenames. Required if SNPs are to be considered and a snpfile is not supplied.')
    parser.add_argument('--gff', type = str, help = 'Genome annotation in gff format.')
    parser.add_argument('--gfftype', type = str, help = 'Source of genome annotation file.', choices = ['GENCODE', 'Ensembl'], required = True)
    parser.add_argument('--genomeFasta', type = str, help = 'Genome sequence in fasta format. Required if SNPs are to be considered.')
    parser.add_argument('--nproc', type = int, help = 'Number of processors to use. Default is 1.', default = 1)
    parser.add_argument('--useSNPs', action = 'store_true', help = 'Consider SNPs?')
    parser.add_argument('--snpfile', type = str, help = 'VCF file of snps to mask. If --useSNPs but a --snpfile is not supplied, a VCF of snps will be created using --controlsamples.')
    parser.add_argument('--maskbed', help = 'Optional. Bed file of positions to mask from analysis.', default = None)
    parser.add_argument('--ROIbed', help = 'Optional. Bed file of specific regions of interest in which to quantify conversions. If supplied, only conversions in these regions will be quantified.', default = None)
    parser.add_argument('--SNPcoverage', type = int, help = 'Minimum coverage to call SNPs. Default = 20', default = 20)
    parser.add_argument('--SNPfreq', type = float, help = 'Minimum variant frequency to call SNPs. Default = 0.4', default = 0.4)
    parser.add_argument('--onlyConsiderOverlap', action = 'store_true', help = 'Only consider conversions seen in both reads of a read pair? Only possible with paired end data.')
    parser.add_argument('--use_g_t', action = 'store_true', help = 'Consider G->T conversions?')
    parser.add_argument('--use_g_c', action = 'store_true', help = 'Consider G->C conversions?')
    parser.add_argument('--use_g_x', action='store_true', help='Consider G->deletion conversions?')
    parser.add_argument('--use_ng_xg', action='store_true', help='Consider NG->deletionG conversions?')
    parser.add_argument('--use_read1', action = 'store_true', help = 'Use read1 when looking for conversions? Only useful with paired end data.')
    parser.add_argument('--use_read2', action = 'store_true', help = 'Use read2 when looking for conversions? Only useful with paired end data.')
    parser.add_argument('--minMappingQual', type = int, help = 'Minimum mapping quality for a read to be considered in conversion counting. STAR unique mappers have MAPQ 255.', required = True)
    parser.add_argument('--minPhred', type = int, help = 'Minimum phred quality score for a base to be considered. Default = 30', default = 30)
    parser.add_argument('--nConv', type = int, help = 'Minimum number of required G->T and/or G->C conversions in a read pair in order for those conversions to be counted. Default is 1.', default = 1)
    parser.add_argument('--outputDir', type = str, help = 'Output directory.', required = True)
    args = parser.parse_args()

    #What type of gff are we working with?
    if args.gfftype == 'GENCODE':
        from assignreads_salmon import getpostmasterassignments, assigntotxs, collapsetogene, readspergene, writeOutput
    elif args.gfftype == 'Ensembl':
        from assignreads_salmon_ensembl import getpostmasterassignments, assigntotxs, collapsetogene, readspergene, writeOutput

    #If we have single end data, considering overlap of paired reads or only one read doesn't make sense
    if args.datatype == 'single':
        args.onlyConsiderOverlap = False
        args.use_read1 = False
        args.use_read2 = False
    
    #Store command line arguments
    suppliedargs = {}
    for arg in vars(args):
        if arg != 'samplenames':
            suppliedargs[arg] = getattr(args, arg)

    #Take in list of samplenames to run pigpen on
    #Derive quant.sf, STAR bams, and postmaster bams
    samplenames = args.samplenames.split(',')
    salmonquants = [os.path.join(x, 'salmon', '{0}.quant.sf'.format(x)) for x in samplenames]
    starbams = [os.path.join(x, 'STAR', '{0}Aligned.sortedByCoord.out.bam'.format(x)) for x in samplenames] #non-deduplicated bams
    #starbams = [os.path.join(x, 'STAR', '{0}.dedup.bam'.format(x)) for x in samplenames]
    postmasterbams = [os.path.join(x, 'postmaster', '{0}.postmaster.bam'.format(x)) for x in samplenames]

    #Take in list of control samples, make list of their corresponding star bams for SNP calling
    if args.controlsamples:
        controlsamples = args.controlsamples.split(',')
        controlindicies = []
    samplebams = []
    controlsamplebams = []
    for ind, x in enumerate(samplenames):
        samplebams.append(starbams[ind])
        if args.controlsamples and x in controlsamples:
            controlsamplebams.append(starbams[ind])

    #We have to be either looking for G->T or G->C, if not both
    if not args.use_g_t and not args.use_g_c:
        print('We have to either be looking for G->T or G->C, if not both! Add argument --use_g_t and/or --use_g_c.')
        sys.exit()

    #We have to be using either read1 or read2 if not both
    if not args.use_read1 and not args.use_read2 and args.datatype == 'paired':
        print('We need to use read1 or read2, if not both! Add argument --use_read1 and/or --use_read2.')
        sys.exit()

    #If we want to only consider overlap, we have to be using both read1 and read2
    if args.onlyConsiderOverlap and (not args.use_read1 or not args.use_read2):
        print('If we are only going to consider overlap between paired reads, we must use both read1 and read2.')
        sys.exit()

    #Make vcf file for snps
    if args.snpfile:
        snps = recordSNPs(args.snpfile)
    if args.useSNPs and not args.snpfile and not args.controlsamples:
        print('ERROR: If we want to consider snps we either have to give control samples for finding snps or a vcf file of snps we already know!')
        sys.exit()
    if args.useSNPs and not args.snpfile:
        if not os.path.exists('snps'):
            os.mkdir('snps')
        vcfFileNames = getSNPs(controlsamplebams, args.genomeFasta, args.SNPcoverage, args.SNPfreq)
        for f in vcfFileNames:
            csi = f + '.csi'
            log = f[:-3] + '.log'
            #Move files to snps directory
            os.rename(f, os.path.join('snps', f))
            os.rename(csi, os.path.join('snps', csi))
            os.rename(log, os.path.join('snps', log))

        os.rename('merged.vcf', os.path.join('snps', 'merged.vcf'))
        os.rename('vcfconcat.log', os.path.join('snps', 'vcfconcat.log'))
        snps = recordSNPs(os.path.join('snps', 'merged.vcf'))

    elif not args.useSNPs and not args.snpfile:
        snps = None

    #Get positions to manually mask if given
    if args.maskbed:
        print('Getting positions to manually mask...')
        maskpositions = readmaskbed(args.maskbed)
    elif not args.maskbed:
        maskpositions = None

    #If there is no supplied bedfile of regions of interest,
    #for each sample, identify conversions, assign conversions to transcripts,
    #and collapse transcript-level measurements to gene-level measurements.
    if not args.ROIbed:
        for ind, sample in enumerate(samplenames):
            #Create parameter dictionary that is unique to this sample
            sampleparams = suppliedargs
            sampleparams['sample'] = sample

            print('Running PIGPEN for {0}...'.format(sample))

            samplebam = samplebams[ind]
            sampleparams['samplebam'] = os.path.abspath(samplebam)
            if args.nproc == 1:
                if args.datatype == 'paired':
                    convs, readcounter = iteratereads_pairedend(samplebam, args.onlyConsiderOverlap, args.use_g_t, args.use_g_c, args.use_g_x, args.use_ng_xg,
                                                                args.use_read1, args.use_read2, args.nConv, args.minMappingQual, args.minPhred, snps, maskpositions, 'high')
                elif args.datatype == 'single':
                    convs, readcounter = iterratereads_singleend(
                        samplebam, args.use_g_t, args.use_g_c, args.nConv, args.minMappingQual, snps, maskpositions, 'high')
            elif args.nproc > 1:
                convs = getmismatches(args.datatype, samplebam, args.onlyConsiderOverlap, snps, maskpositions, args.nConv,
                                      args.minMappingQual, args.nproc, args.use_g_t, args.use_g_c, args.use_g_x, args.use_ng_xg, args.use_read1, args.use_read2, args.minPhred)

            print('Getting posterior probabilities from salmon alignment file...')
            postmasterbam = postmasterbams[ind]
            sampleparams['postmasterbam'] = os.path.abspath(postmasterbam)
            pprobs = getpostmasterassignments(postmasterbam)
            print('Assinging conversions to transcripts...')
            txconvs = assigntotxs(pprobs, convs)
            print('Collapsing transcript level conversion counts to gene level...')
            tx2gene, geneid2genename, geneconvs = collapsetogene(txconvs, args.gff)
            print('Counting number of reads assigned to each gene...')
            salmonquant = salmonquants[ind]
            sampleparams['salmonquant'] = os.path.abspath(salmonquant)
            genecounts = readspergene(salmonquant, tx2gene)
            print('Writing output...')
            if not os.path.exists(args.outputDir):
                os.mkdir(args.outputDir)
            outputfile = os.path.join(args.outputDir, sample + '.pigpen.txt')
            writeOutput(sampleparams, geneconvs, genecounts, geneid2genename, outputfile, args.use_g_t, args.use_g_c, args.use_g_x, args.use_ng_xg)
            print('Done!')

    #If there is a bed file of regions of interest supplied, then use that. Don't use the salmon/postmaster quantifications.
    elif args.ROIbed:
        #Make fasta index
        command = ['samtools', 'faidx', args.genomeFasta]
        subprocess.call(command)
        faidx = args.genomeFasta + '.fai'

        #Create chrsort
        command = ['cut', '-f' '1,2', faidx]
        with open('chrsort.txt', 'w') as outfh:
            subprocess.run(command, stdout = outfh)

        for ind, sample in enumerate(samplenames):
            #Create parameter dictionary that is unique to this sample
            sampleparams = suppliedargs
            sampleparams['sample'] = sample

            print('Running PIGPEN for {0}...'.format(sample))

            samplebam = samplebams[ind]
            sampleparams['samplebam'] = os.path.abspath(samplebam)
            if args.nproc == 1:
                if args.datatype == 'paired':
                    convs, readcounter = iteratereads_pairedend(samplebam, args.onlyConsiderOverlap, args.use_g_t, args.use_g_c, args.use_g_x, args.use_ng_xg,
                                                                args.use_read1, args.use_read2, args.nConv, args.minMappingQual, args.minPhred, snps, maskpositions, 'high')
                elif args.datatype == 'single':
                    convs, readcounter = iterratereads_singleend(
                        samplebam, args.use_g_t, args.use_g_c, args.nConv, args.minMappingQual, snps, maskpositions, 'high')

            elif args.nproc > 1:
                convs = getmismatches(args.datatype, samplebam, args.onlyConsiderOverlap, snps, maskpositions,
                                      args.nConv, args.minMappingQual, args.nproc, args.use_g_t, args.use_g_c, args.use_g_x, args.use_ng_xg, args.use_read1, args.use_read2, args.minPhred)

            print('Assigning reads to genes in supplied bed file...')
            overlaps, numpairs = getReadOverlaps(samplebam, args.ROIbed, 'chrsort.txt')
            read2gene = processOverlaps(overlaps, numpairs)
            numreadspergene, convsPerGene = getPerGene(convs, read2gene)
            if not os.path.exists(args.outputDir):
                os.mkdir(args.outputDir)
            outputfile = os.path.join(args.outputDir, sample + '.pigpen.txt')
            writeConvsPerGene(sampleparams, numreadspergene, convsPerGene, outputfile, args.use_g_t, args.use_g_c, args.use_g_x, args.use_ng_xg)