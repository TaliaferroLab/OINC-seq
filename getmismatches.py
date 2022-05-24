import pysam
import os
import re
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
from snps import getSNPs, recordSNPs
import pickle
import multiprocessing as mp
import subprocess


def revcomp(nt):
    revcompdict = {
        'G' : 'C',
        'C' : 'G',
        'A' : 'T',
        'T' : 'A',
        'N' : 'N',
        'g' : 'c',
        'c' : 'g',
        'a' : 't',
        't' : 'a',
        'n' : 'n'
    }

    nt_rc = revcompdict[nt]

    return nt_rc


def iteratereads_singleend(bam, snps = None):
    #Read through a bam containing single end reads (or if it contains paired end reads, just use read 1)
    #Find nt conversion locations for each read.
    #Store the number of each conversion for each read in a dictionary.

    #Bam must contain MD tag.

    readcounter = 0
    convs = {} #{readid : dictionary of all conversions}
    with pysam.AlignmentFile(bam, 'r') as infh:
        print('Finding nucleotide conversions in {0}...'.format(os.path.basename(bam)))
        for read in infh.fetch(until_eof = True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            if read.is_read1:
                readcounter +=1
                if readcounter % 10000000 == 0:
                    print('Finding nucleotide conversions in read {0}...'.format(readcounter))

                #Check mapping quality
                #For nextgenmap, max mapq is 60
                if read.mapping_quality < 255:
                    continue
                
                queryname = read.query_name
                queryseq = read.query_sequence #this is always on the + strand, no matter what strand the read maps to
                chrom = read.reference_name
                qualities = list(read.query_qualities)

                #Get a set of snp locations if we have them
                if snps:
                    if chrm in snps:
                        snplocations = snps[chrm] #set of coordinates to mask
                    else:
                        snplocations = None
                else:
                    snplocations = None

                if read.is_reverse:
                    strand = '-'
                elif not read.is_reverse:
                    strand = '+'

                alignedpairs = read.get_aligned_pairs(with_seq = True)
                readqualities = list(read.query_qualities)
                convs_in_read = getmismatches_singleend(alignedpairs, queryseq, readqualities, strand, chrom, snplocations)
                
                convs[queryname] = convs_in_read


    #Pickle and write convs?
    return convs
                

def getmismatches_singleend(alignedpairs, queryseq, readqualities, strand, chrom, snplocations):
    #remove tuples that have None
    #These are either intronic or might have been soft-clipped
    #Tuples are (querypos, (0-based) refpos, refsequence)
    #If there is a substitution, refsequence is lower case

    #refnt and querynt, as supplied by pysam, are always + strand
    #everything here is assuming read is on sense strand

    #snplocations is a set of chrm_coord locations. At these locations, all queries will be treated as not having a conversion

    alignedpairs = [x for x in alignedpairs if None not in x]

    #if we have snps, remove their locations from alignedpairs
    #snplocations is a set of 0-based coordinates of snp locations to mask
    if snplocations:
        alignedpairs = [x for x in alignedpairs if x[1] not in snplocations]
    
    convs = {} #counts of conversions x_y where x is reference sequence and y is query sequence
    convlocations = defaultdict(list) #{type of conv : [locations of conversion]}

    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n']

    #initialize dictionary
    for conv in possibleconvs:
        convs[conv] = 0

    for alignedpair in alignedpairs:
        refnt = alignedpair[2]
        
        #if reference is N, skip this position
        if refnt == 'N' or refnt == 'n':
            continue

        if strand == '-':
            refnt = revcomp(refnt)

        #Not a conversion
        if refnt.isupper():
            conv = refnt.lower() + '_' + refnt.lower()
        #Is a conversion
        elif refnt.islower():
            querynt = queryseq[alignedpair[0]]
            if strand == '-':
                querynt = revcomp(querynt)
            conv = refnt.lower() + '_' + querynt.lower()

        #If the quality at this position passes threshold, record the conversion.
        #Otherwise, skip it.
        qscore = readqualities[alignedpair[0]]
        if qscore >= 30: #can change this later
            convs[conv] +=1
        else:
            pass

    return convs

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.mate_is_unmapped or read.is_unmapped:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def findsnps(controlbams, genomefasta, minCoverage = 20, minVarFreq = 0.02):
    controlbams = controlbams.split(',')
    getSNPs(controlbams, genomefasta, minCoverage, minVarFreq)
    snps = recordSNPs('merged.vcf')

    return snps


def iteratereads_pairedend(bam, onlyConsiderOverlap, use_g_t, use_g_c, nConv, snps=None, maskpositions=None, verbosity='high'):
    #Iterate over reads in a paired end alignment file.
    #Find nt conversion locations for each read.
    #For locations interrogated by both mates of read pair, conversion must exist in both mates in order to count
    #Store the number of conversions for each read in a dictionary

    #Quality score array is always in the same order as query_sequence, which is always on the + strand
    #Bam must contain MD tags

    if onlyConsiderOverlap == 'True':
        onlyConsiderOverlap = True
    elif onlyConsiderOverlap == 'False':
        onlyConsiderOverlap = False
        
    queriednts = []
    readcounter = 0
    convs = {} #{readid : dictionary of all conversions}
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'r') as infh:
        if verbosity == 'high':
            print('Finding nucleotide conversions in {0}...'.format(os.path.basename(bam)))
        for read1, read2 in read_pair_generator(infh):
            
            #Just double check that the pairs are matched
            if read1.query_name != read2.query_name:
                continue

            #Check mapping quality
            #MapQ is 255 for uniquely aligned reads FOR STAR ONLY
            if read1.mapping_quality < 255 or read2.mapping_quality < 255:
                continue

            readcounter +=1
            if readcounter % 1000000 == 0:
                if verbosity == 'high':
                    print('Finding nucleotide conversions in read {0}...'.format(readcounter))
                
            queryname = read1.query_name
            chrm = read1.reference_name

            #Get a set of positions to mask (snps + locations we want to mask)
            #Get a set of snp locations if we have them
            if snps:
                if chrm in snps:
                    snplocations = snps[chrm] #set of snp coordinates to mask
                else:
                    snplocations = None
            else:
                snplocations = None

            #Get a set of locations to mask if we have them
            if maskpositions:
                if chrm in maskpositions:
                    masklocations = maskpositions[chrm] #set of coordinates to manually mask
                else:
                    masklocations = None
            else:
                masklocations = None

            #combine snps and manually masked positions into one set
            #this combined set will be masklocations
            if snplocations and masklocations:
                masklocations.update(snplocations)
            elif snplocations and not masklocations:
                masklocations = snplocations
            elif masklocations and not snplocations:
                masklocations = masklocations

            read1queryseq = read1.query_sequence
            read1alignedpairs = read1.get_aligned_pairs(with_seq = True)
            if read1.is_reverse:
                read1strand = '-'
            elif not read1.is_reverse:
                read1strand = '+'

            read2queryseq = read2.query_sequence
            read2alignedpairs = read2.get_aligned_pairs(with_seq = True)
            if read2.is_reverse:
                read2strand = '-'
            elif not read2.is_reverse:
                read2strand = '+'

            read1qualities = list(read1.query_qualities) #phred scores
            read2qualities = list(read2.query_qualities)

            convs_in_read = getmismatches_pairedend(read1alignedpairs, read2alignedpairs, read1queryseq, read2queryseq, read1qualities, read2qualities, read1strand, read2strand, masklocations, onlyConsiderOverlap, nConv, use_g_t, use_g_c)
            queriednts.append(sum(convs_in_read.values()))
            convs[queryname] = convs_in_read

    if verbosity == 'high':
        print('Queried {0} read pairs in {1}.'.format(readcounter, os.path.basename(bam)))

    pysam.set_verbosity(save)
    #Pickle and write convs?
    return convs, readcounter


def getmismatches_pairedend(read1alignedpairs, read2alignedpairs, read1queryseq, read2queryseq, read1qualities, read2qualities, read1strand, read2strand, masklocations, onlyoverlap, nConv, use_g_t, use_g_c):
    #remove tuples that have None
    #These are either intronic or might have been soft-clipped
    #Tuples are (querypos, refpos, refsequence)
    #If there is a substitution, refsequence is lower case

    #remove positions where querypos is None
    #i'm pretty sure these query positions won't have quality scores
    read1alignedpairs = [x for x in read1alignedpairs if x[0] != None]
    read2alignedpairs = [x for x in read2alignedpairs if x[0] != None]

    #Add quality scores to alignedpairs tuples
    #will now be (querypos, refpos, refsequence, qualityscore)
    read1ap_withq = []
    for ind, x in enumerate(read1alignedpairs):
        x += (read1qualities[ind],)
        read1ap_withq.append(x)
    read1alignedpairs = read1ap_withq
    
    read2ap_withq = []
    for ind, x in enumerate(read2alignedpairs):
        x += (read2qualities[ind],)
        read2ap_withq.append(x)
    read2alignedpairs = read1ap_withq

    #Now remove positions where refsequence is None
    #These may be places that got soft-clipped
    read1alignedpairs = [x for x in read1alignedpairs if None not in x]
    read2alignedpairs = [x for x in read2alignedpairs if None not in x]

    #if we have locations to mask, remove their locations from read1alignedpairs and read2alignedpairs
    #masklocations is a set of 0-based coordinates of snp locations to mask
    if masklocations:
        read1alignedpairs = [x for x in read1alignedpairs if x[1] not in masklocations]
        read2alignedpairs = [x for x in read2alignedpairs if x[1] not in masklocations]
    
    convs = {} #counts of conversions x_y where x is reference sequence and y is query sequence

    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n']

    #initialize dictionary
    for conv in possibleconvs:
        convs[conv] = 0

    #For locations interrogated by both mates of read pair, conversion must exist in both mates in order to count
    #These locations (as defined by their reference positions) would be found both in read1alignedpairs and read2alignedpairs

    #Get the ref positions queried by the two reads
    r1dict = {} #{reference position : [queryposition, reference sequence, quality]}
    r2dict = {}
    for x in read1alignedpairs:
        r1dict[int(x[1])] = [x[0], x[2], x[3]]
    for x in read2alignedpairs:
        r2dict[int(x[1])] = [x[0], x[2], x[3]]

    mergedalignedpairs = {} # {refpos : [R1querypos, R2querypos, R1refsequence, R2refsequence, R1quality, R2quality]}
    #For positions only in R1 or R2, querypos and refsequence are NA for the other read
    for refpos in r1dict:
        r1querypos = r1dict[refpos][0]
        r1refseq = r1dict[refpos][1]
        r1quality = r1dict[refpos][2]
        if refpos in mergedalignedpairs: #this should not be possible because we are looking at r1 first
            r2querypos = mergedalignedpairs[refpos][1]
            r2refseq = mergedalignedpairs[refpos][3]
            r2quality = mergedalignedpairs[refpos][5]
            mergedalignedpairs[refpos] = [r1querypos, r2querypos, r1refseq, r2refseq, r1quality, r2quality]
        else:
            mergedalignedpairs[refpos] = [r1querypos, 'NA', r1refseq, 'NA', r1quality, 'NA']

    for refpos in r2dict:
        #same thing
        r2querypos = r2dict[refpos][0]
        r2refseq = r2dict[refpos][1]
        r2quality = r2dict[refpos][2]
        if refpos in mergedalignedpairs: #if we saw it for r1
            r1querypos = mergedalignedpairs[refpos][0]
            r1refseq = mergedalignedpairs[refpos][2]
            r1quality = mergedalignedpairs[refpos][4]
            mergedalignedpairs[refpos] = [r1querypos, r2querypos, r1refseq, r2refseq, r1quality, r2quality]
        else:
            mergedalignedpairs[refpos] = ['NA', r2querypos, 'NA', r2refseq, 'NA', r2quality]


    #Now go through mergedalignedpairs, looking for conversions.
    #For positions observed both in r1 and r2, a conversion must be present in both reads,
    #otherwise it will be recorded as not having a conversion.

    for refpos in mergedalignedpairs:
        r1querypos = mergedalignedpairs[refpos][0]
        r2querypos = mergedalignedpairs[refpos][1]
        r1refseq = mergedalignedpairs[refpos][2]
        r2refseq = mergedalignedpairs[refpos][3]
        r1quality = mergedalignedpairs[refpos][4]
        r2quality = mergedalignedpairs[refpos][5]

        if r1querypos != 'NA' and r2querypos == 'NA': #this position queried by r1 only
            if read1strand == '-':
                #refseq needs to equal the sense strand (it is always initially defined as the + strand). read1 is always the sense strand.
                r1refseq = revcomp(r1refseq)

            #If reference is N, skip this position
            if r1refseq == 'N' or r1refseq == 'n':
                continue

            if r1refseq.isupper(): #not a conversion
                conv = r1refseq.lower() + '_' + r1refseq.lower()
            elif r1refseq.islower(): #is a conversion
                querynt = read1queryseq[r1querypos]
                if read1strand == '-':
                    querynt = revcomp(querynt)
                conv = r1refseq.lower() + '_' + querynt.lower()

            if r1quality >= 30 and onlyoverlap == False:
                convs[conv] +=1
            else:
                pass

        elif r1querypos == 'NA' and r2querypos != 'NA': #this position is queried by r2 only
            if read1strand == '-':
                r2refseq = revcomp(r2refseq) #reference seq is independent of which read we are talking about

            if r2refseq == 'N' or r2refseq == 'n':
                continue

            if r2refseq.isupper(): #not a conversion
                conv = r2refseq.lower() + '_' + r2refseq.lower()
            elif r2refseq.islower(): #is a conversion
                querynt = read2queryseq[r2querypos]
                if read1strand == '-':
                    #Read1 is always the sense strand. r1queryseq and r2queryseq are always + strand
                    #The reference sequence is always on the + strand.
                    #If read1 is on the - strand, we have already flipped reference seq (see a few lines above).
                    #We need to flip read2queryseq so that it is also - strand.
                    querynt = revcomp(querynt)
                conv = r2refseq.lower() + '_' + querynt.lower()

            if r2quality >= 30 and onlyoverlap == False:
                convs[conv] +=1
            else:
                pass

        elif r1querypos != 'NA' and r2querypos != 'NA': #this position is queried by both reads
            if read1strand == '-':
                r1refseq = revcomp(r1refseq)
                r2refseq = revcomp(r2refseq)
            
            if r1refseq == 'N' or r2refseq == 'N' or r1refseq == 'n' or r2refseq == 'n':
                continue

            #If the position is not high quality in both r1 and r2, skip it
            if r1quality < 30 and r2quality < 30:
                continue

            if r1refseq.isupper() and r2refseq.isupper(): #both reads agree it is not a conversion
                conv = r1refseq.lower() + '_' + r1refseq.lower()
                convs[conv] += 1

            elif r1refseq.isupper() and r2refseq.islower(): #r1 says no conversion, r2 says conversion, so we say no conversion
                conv = r1refseq.lower() + '_' + r1refseq.lower()
                convs[conv] += 1

            elif r1refseq.islower() and r2refseq.isupper(): #r1 says conversion, r2 says no conversion, so we say no conversion
                conv = r2refseq.lower() + '_' + r2refseq.lower()
                convs[conv] += 1

            elif r1refseq.islower() and r2refseq.islower(): #both reads say there was a conversion
                r1querynt = read1queryseq[r1querypos]
                r2querynt = read2queryseq[r2querypos]
                if read1strand == '-':
                    r1querynt = revcomp(r1querynt)
                    r2querynt = revcomp(r2querynt)

                #If the query nts don't match, skip this position
                if r1querynt == r2querynt:
                    conv = r1refseq.lower() + '_' + r1querynt.lower()
                    convs[conv] +=1
                else:
                    pass

    #Does the number of g_t and/or g_c conversions meet our threshold?
    if use_g_t and use_g_c:
        if convs['g_t'] + convs['g_c'] >= nConv:
            pass
        elif convs['g_t'] + convs['g_c'] < nConv:
            convs['g_t'] = 0
            convs['g_c'] = 0
    elif use_g_t and not use_g_c:
        if convs['g_t'] >= nConv:
            pass
        elif convs['g_t'] < nConv:
            convs['g_t'] = 0
            convs['g_c'] = 0
    elif use_g_c and not use_g_t:
        if convs['g_c'] >= nConv:
            pass
        elif convs['g_c'] < nConv:
            convs['g_c'] = 0
            convs['g_t'] = 0
    elif not use_g_t and not use_g_c:
        print('ERROR: we have to be looking for at least either G->T or G->C if not both!!')
        sys.exit()

    return convs


def summarize_convs(convs, outfile):
    #Take in a dictionary of conversions and calculate conversion rates for all conversions
    #convs has the form {readid : {a_a : count, a_t : count, etc.}}

    a = 0 #running counter of all *reference* a's (converted and not converted)
    a_t = 0 #running counter of all a to t conversions
    a_c = 0
    a_g = 0
    a_n = 0
    c = 0
    c_t = 0
    c_g = 0
    c_a = 0
    c_n = 0
    g = 0
    g_t = 0
    g_a = 0
    g_c = 0
    g_n = 0
    t = 0
    t_a = 0
    t_g = 0
    t_c = 0
    t_n = 0

    for read in convs:
        conv_in_read = convs[read]
        a += (conv_in_read['a_a'] + conv_in_read['a_t'] + conv_in_read['a_c'] + conv_in_read['a_g'] + conv_in_read['a_n'])
        a_t += conv_in_read['a_t']
        a_c += conv_in_read['a_c']
        a_g += conv_in_read['a_g']
        a_n += conv_in_read['a_n']

        c += (conv_in_read['c_a'] + conv_in_read['c_t'] + conv_in_read['c_c'] + conv_in_read['c_g'] + conv_in_read['c_n'])
        c_t += conv_in_read['c_t']
        c_a += conv_in_read['c_a']
        c_g += conv_in_read['c_g']
        c_n += conv_in_read['c_n']

        g += (conv_in_read['g_a'] + conv_in_read['g_t'] + conv_in_read['g_c'] + conv_in_read['g_g'] + conv_in_read['g_n'])
        g_t += conv_in_read['g_t']
        g_a += conv_in_read['g_a']
        g_c += conv_in_read['g_c']
        g_n += conv_in_read['g_n']

        t += (conv_in_read['t_a'] + conv_in_read['t_t'] + conv_in_read['t_c'] + conv_in_read['t_g'] + conv_in_read['t_n'])
        t_g += conv_in_read['t_g']
        t_a += conv_in_read['t_a']
        t_c += conv_in_read['t_c']
        t_n += conv_in_read['t_n']

    totalnt = a + c + g + t
    totalconv = a_t + a_c + a_g + c_t + c_a + c_g + g_t + g_a + g_c + t_g + t_a + t_c

    totalerrorrate = totalconv / totalnt

    with open(outfile, 'w') as outfh:
        outfh.write(('\t').join([
            'Acount', 'A_Tcount', 'A_Ccount', 'A_Gcount', 'A_Ncount',
            'Ccount', 'C_Tcount', 'C_Acount', 'C_Gcount', 'C_Ncount',
            'Gcount', 'G_Tcount', 'G_Ccount', 'G_Acount', 'G_Ncount',
            'Tcount', 'T_Acount', 'T_Ccount', 'T_Gcount', 'T_Ncount',
            'A_Trate', 'A_Crate', 'A_Grate', 'A_Nrate',
            'C_Trate', 'C_Arate', 'C_Grate', 'C_Nrate',
            'G_Trate', 'G_Crate', 'G_Arate', 'G_Nrate',
            'T_Arate', 'T_Crate', 'T_Grate', 'T_Nrate', 'totalnt', 'totalconv', 'totalerrorrate'
        ]) + '\n')

        outfh.write(('\t').join([
            str(a), str(a_t), str(a_c), str(a_g), str(a_n),
            str(c), str(c_t), str(c_a), str(c_g), str(c_n),
            str(g), str(g_t), str(g_c), str(g_a), str(g_n),
            str(t), str(t_a), str(t_c), str(t_g), str(t_n),
            str(a_t / a), str(a_c / a), str(a_g / a), str(a_n / a),
            str(c_t / c), str(c_a / c), str(c_g / c), str(c_n / c),
            str(g_t / g), str(g_c / g), str(g_a / g), str(g_n / g),
            str(t_a / t), str(t_c / t), str(t_g / t), str(t_n / t),
            str(totalnt), str(totalconv), str(totalerrorrate)
        ]))

def split_bam(bam, nproc):
    #Split bam using samtools instead of bamtools (may be faster)
    #Check for index
    if os.path.exists(os.path.abspath(bam) + '.bai'):
        pass
    else:
        indexCMD = 'samtools index ' + os.path.abspath(bam)
        index = subprocess.Popen(indexCMD, shell = True)
        index.wait()

    #Get chromosome names
    idxstatsCMD = 'samtools idxstats ' + os.path.abspath(bam)
    idxstats = subprocess.Popen(idxstatsCMD, shell = True, stdout = subprocess.PIPE)
    chrms = []
    for line in idxstats.stdout:
        line = line.decode('UTF-8').strip().split('\t')
        chrm = line[0]
        if chrm != '*':
            chrms.append(chrm)

    splitbams = []
    for chrm in chrms:
        filebasename = chrm + '_SPLIT.bam'
        filepath = os.path.join(os.path.dirname(os.path.abspath(bam)), filebasename)
        splitbams.append(filepath)
        splitCMD = 'samtools view -@ ' + str(nproc) + ' -b ' + os.path.abspath(bam) + ' ' + chrm + ' > ' + filepath
        s = subprocess.Popen(splitCMD, shell = True)
        s.wait()

    return splitbams


def getmismatches(bam, onlyConsiderOverlap, snps, maskpositions, nConv, nproc, use_g_t, use_g_c):
    #Actually run the mismatch code (calling iteratereads_pairedend)
    #use multiprocessing
    #If there's only one processor, easier to use iteratereads_pairedend() directly.

    pool = mp.Pool(processes = int(nproc))
    print('Using {0} processors to identify mismatches in {1}.'.format(nproc, bam))
    splitbams = split_bam(bam, int(nproc))
    argslist = []
    for x in splitbams:
        argslist.append((x, bool(onlyConsiderOverlap), bool(
            use_g_t), bool(use_g_c), nConv, snps, maskpositions, 'low'))

    #items returned from iteratereads_pairedend are in a list, one per process
    totalreadcounter = 0 #number of reads across all the split bams
    results = pool.starmap(iteratereads_pairedend, argslist) #thhis actually returns two things, convs and readcounter
    #so i bet this is a nested list where the first item in each list in a convs and the second item is a readcounter
    convs_split = []
    for result in results:
        convs_split.append(result[0])
    for result in results:
        totalreadcounter += result[1]

    print('Queried {0} read pairs in {1}.'.format(totalreadcounter, os.path.basename(bam)))

    #Reorganize convs_split into convs as it is without multiprocessing
    convs = {} #{readid : dictionary of all conversions}
    for splitconv in convs_split:
        convs.update(splitconv)

    #cleanup
    for splitbam in splitbams:
        os.remove(splitbam)

    return convs


    

        
if __name__ == '__main__':
    iteratereads_pairedend(sys.argv[1], False, True, True, 1, None, 'high')

    
