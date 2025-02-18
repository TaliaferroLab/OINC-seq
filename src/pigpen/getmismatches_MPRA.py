import pysam
import os
import re
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import subprocess
import argparse


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
        'n' : 'n', 
        'X': 'X',
        None: None,
        'NA': 'NA'
    }

    nt_rc = revcompdict[nt]

    return nt_rc


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


def iteratereads_pairedend(bam, onlyConsiderOverlap, use_read1, use_read2, nConv, minMappingQual, minPhred=30):
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
    convs = {}  # {oligoname : [readcount, {dictionary of all conversions}]
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'r') as infh:
        print('Finding nucleotide conversions in {0}...'.format(
            os.path.basename(bam)))
        for read1, read2 in read_pair_generator(infh):

            #Just double check that the pairs are matched
            if read1.query_name != read2.query_name:
                continue

            if read1.reference_name != read2.reference_name:
                continue

            #Check mapping quality
            #MapQ is 255 for uniquely aligned reads FOR STAR ONLY
            #MapQ is (I think) >= 2 for uniquely aligned reads from bowtie2
            if read1.mapping_quality < minMappingQual or read2.mapping_quality < minMappingQual:
                continue

            readcounter += 1
            if readcounter % 1000000 == 0:
                print('Finding nucleotide conversions in read {0}...'.format(
                    readcounter))

            oligo = read1.reference_name
            read1queryseq = read1.query_sequence
            read1alignedpairs = read1.get_aligned_pairs(with_seq=True)
            if read1.is_reverse:
                read1strand = '-'
            elif not read1.is_reverse:
                read1strand = '+'

            read2queryseq = read2.query_sequence
            read2alignedpairs = read2.get_aligned_pairs(with_seq=True)
            if read2.is_reverse:
                read2strand = '-'
            elif not read2.is_reverse:
                read2strand = '+'

            read1qualities = list(read1.query_qualities)  # phred scores
            read2qualities = list(read2.query_qualities)

            convs_in_read = getmismatches_pairedend(read1alignedpairs, read2alignedpairs, read1queryseq, read2queryseq, read1qualities, read2qualities, read1strand, read2strand, onlyConsiderOverlap, nConv, minPhred, use_read1, use_read2)
            queriednts.append(sum(convs_in_read.values()))
            
            #Add to our running dictionary
            if oligo not in convs:
                convs[oligo] = [1, convs_in_read]
            elif oligo in convs:
                readcount = convs[oligo][0]
                readcount +=1
                oligodict = convs[oligo][1]
                for conv in convs_in_read:
                    convcount = convs_in_read[conv]
                    oligodict[conv] += convcount
                convs[oligo] = [readcount, oligodict]

    return convs


def getmismatches_pairedend(read1alignedpairs, read2alignedpairs, read1queryseq, read2queryseq, read1qualities, read2qualities, read1strand, read2strand, masklocations, onlyoverlap, nConv, minPhred, use_read1, use_read2):
    #remove tuples that have None
    #These are either intronic or might have been soft-clipped
    #Tuples are (querypos, refpos, refsequence)
    #If there is a substitution, refsequence is lower case

    #For now, forget insertions. Get rid of any position where reference position is None.
    read1alignedpairs = [x for x in read1alignedpairs if x[1] != None]
    read2alignedpairs = [x for x in read2alignedpairs if x[1] != None]
    read1alignedpairs = [x for x in read1alignedpairs if x[2] != None]
    read2alignedpairs = [x for x in read2alignedpairs if x[2] != None]

    #Add quality scores and query sequences
    #will now be (querypos, refpos, refsequence, querysequence, qualscore)
    for x in range(len(read1alignedpairs)):
        alignedpair = read1alignedpairs[x]
        querypos = alignedpair[0]
        if querypos != None:
            querynt = read1queryseq[querypos]
            qualscore = read1qualities[querypos]
        elif querypos == None:
            querynt = 'X'  # there is no query nt for a deletion
            qualscore = 37  # there's no query position here, so make up a quality score
        alignedpair = alignedpair + (querynt, qualscore)
        read1alignedpairs[x] = alignedpair

    for x in range(len(read2alignedpairs)):
        alignedpair = read2alignedpairs[x]
        querypos = alignedpair[0]
        if querypos != None:
            querynt = read2queryseq[querypos]
            qualscore = read2qualities[querypos]
        elif querypos == None:
            querynt = 'X'
            qualscore = 37
        alignedpair = alignedpair + (querynt, qualscore)
        read2alignedpairs[x] = alignedpair

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
        't_a', 't_t', 't_c', 't_g', 't_n',
        'a_x', 'g_x', 'c_x', 't_x', 'ng_xg']

    #initialize dictionary
    for conv in possibleconvs:
        convs[conv] = 0

    #For locations interrogated by both mates of read pair, conversion must exist in both mates in order to count
    #These locations (as defined by their reference positions) would be found both in read1alignedpairs and read2alignedpairs

    #Get the ref positions queried by the two reads
    r1dict = {} #{reference position: [queryposition, reference sequence, querysequence, quality]}
    r2dict = {}
    for x in read1alignedpairs:
        r1dict[int(x[1])] = [x[0], x[2].upper(), x[3].upper(), x[4]]
    for x in read2alignedpairs:
        r2dict[int(x[1])] = [x[0], x[2].upper(), x[3].upper(), x[4]]

    # {refpos : [R1querypos, R2querypos, R1refsequence, R2refsequence, R1querysequence, R2querysequence, R1quality, R2quality]}
    mergedalignedpairs = {}
    #For positions only in R1 or R2, querypos and refsequence are NA for the other read
    for refpos in r1dict:
        r1querypos = r1dict[refpos][0]
        r1refseq = r1dict[refpos][1]
        r1queryseq = r1dict[refpos][2]
        r1quality = r1dict[refpos][3]
        if refpos in mergedalignedpairs: #this should not be possible because we are looking at r1 first
            r2querypos = mergedalignedpairs[refpos][1]
            r2refseq = mergedalignedpairs[refpos][3]
            r2queryseq = mergedalignedpairs[refpos][5]
            r2quality = mergedalignedpairs[refpos][7]
            mergedalignedpairs[refpos] = [r1querypos, r2querypos, r1refseq,
                                          r2refseq, r1queryseq, r2queryseq, r1quality, r2quality]
        else:
            mergedalignedpairs[refpos] = [r1querypos, 'NA', r1refseq, 'NA', r1queryseq, 'NA', r1quality, 'NA']

    for refpos in r2dict:
        #same thing
        r2querypos = r2dict[refpos][0]
        r2refseq = r2dict[refpos][1]
        r2queryseq = r2dict[refpos][2]
        r2quality = r2dict[refpos][3]
        if refpos in mergedalignedpairs: #if we saw it for r1
            r1querypos = mergedalignedpairs[refpos][0]
            r1refseq = mergedalignedpairs[refpos][2]
            r1queryseq = mergedalignedpairs[refpos][4]
            r1quality = mergedalignedpairs[refpos][6]
            mergedalignedpairs[refpos] = [r1querypos, r2querypos, r1refseq,
                                          r2refseq, r1queryseq, r2queryseq, r1quality, r2quality]
        else:
            mergedalignedpairs[refpos] = ['NA', r2querypos, 'NA', r2refseq, 'NA', r2queryseq, 'NA', r2quality]

    #If we are only using read1 or only using read2, replace the positions in the non-used read with NA
    for refpos in mergedalignedpairs:
        r1querypos, r2querypos, r1refseq, r2refseq, r1queryseq, r2queryseq, r1quality, r2quality = mergedalignedpairs[refpos]
        if use_read1 and not use_read2:
            updatedlist = [r1querypos, 'NA', r1refseq,
                           'NA', r1queryseq, 'NA', r1quality, 'NA']
            mergedalignedpairs[refpos] = updatedlist
        elif use_read2 and not use_read1:
            updatedlist = ['NA', r2querypos, 'NA',
                           r2refseq, 'NA', r2queryseq, 'NA', r2quality]
            mergedalignedpairs[refpos] = updatedlist
        elif not use_read1 and not use_read2:
            print('ERROR: we have to use either read1 or read2, if not both.')
            sys.exit()
        elif use_read1 and use_read2:
            pass

    #Now go through mergedalignedpairs, looking for conversions.
    #For positions observed both in r1 and r2, queryseq in both reads must match, otherwise the position is skipped.
    #We are now keeping track of deletions as either g_x (reference G, query deletion) or ng_xg (ref nt 5' of G deleted in query)
    #We have observed that sometimes RT skips the nucleotide *after* an oxidized G (after being from the RT's point of view)

    for refpos in mergedalignedpairs:
        conv = None
        conv2 = None #sometimes we can have 2 convs (for example the first nt of ng_xg could be both g_x and ng_xg)
        r1querypos = mergedalignedpairs[refpos][0]
        r2querypos = mergedalignedpairs[refpos][1]
        r1refseq = mergedalignedpairs[refpos][2]
        r2refseq = mergedalignedpairs[refpos][3]
        r1queryseq = mergedalignedpairs[refpos][4]
        r2queryseq = mergedalignedpairs[refpos][5]
        r1quality = mergedalignedpairs[refpos][6]
        r2quality = mergedalignedpairs[refpos][7]

        if r1querypos != 'NA' and r2querypos == 'NA': #this position queried by r1 only
            if read1strand == '-':
                #refseq needs to equal the sense strand (it is always initially defined as the + strand). read1 is always the sense strand.
                r1refseq = revcomp(r1refseq)
                r1queryseq = revcomp(r1queryseq)

            #If reference is N, skip this position
            if r1refseq == 'N' or r1refseq == 'n':
                continue
            conv = r1refseq.lower() + '_' + r1queryseq.lower()

            if r1queryseq == 'X':
                #Check if there is a reference G downstream of this position
                if read1strand == '+':
                    downstreamrefpos = refpos + 1
                    #It's possible that downstreamrefpos is not in mergedalignedpairs because this position is at the end of the read
                    try:
                        downstreamrefseq = mergedalignedpairs[downstreamrefpos][2].upper()
                        downstreamqueryseq = mergedalignedpairs[downstreamrefpos][4].upper()
                    except KeyError:
                        downstreamrefseq, downstreamqueryseq = None, None
                elif read1strand == '-':
                    downstreamrefpos = refpos - 1
                    try:
                        downstreamrefseq = mergedalignedpairs[downstreamrefpos][2].upper()
                        downstreamqueryseq = mergedalignedpairs[downstreamrefpos][4].upper()
                    except KeyError:
                        downstreamrefseq, downstreamqueryseq = None, None
                    downstreamrefseq = revcomp(downstreamrefseq)
                    downstreamqueryseq = revcomp(downstreamqueryseq)
                if downstreamrefseq == 'G':
                    conv2 = 'ng_xg'
                    #If this is a non-g deletion and is downstream of a g, we can't be sure if this deletion is due to this nucleotide or the downstream g
                    if conv in ['a_x', 't_x', 'c_x']:
                        conv = None

            #Add conv(s) to dictionary
            if r1quality >= minPhred and onlyoverlap == False:
                # there will be some conversions (e.g. a_x that are not in convs)
                if conv in convs:
                    convs[conv] +=1
                if conv2 == 'ng_xg' and conv != 'g_x':
                    convs[conv2] +=1

        elif r1querypos == 'NA' and r2querypos != 'NA': #this position is queried by r2 only
            if read1strand == '-':
                # reference seq is independent of which read we are talking about
                #Read1 is always the sense strand. r1queryseq and r2queryseq are always + strand
                #The reference sequence is always on the + strand.
                #If read1 is on the - strand, we have already flipped reference seq (see a few lines above).
                #We need to flip read2queryseq so that it is also - strand.
                r2refseq = revcomp(r2refseq)
                r2queryseq = revcomp(r2queryseq)

            if r2refseq == 'N' or r2refseq == 'n':
                continue

            conv = r2refseq.lower() + '_' + r2refseq.lower()
            if r2queryseq == 'X':
                #Check if there is a reference G downstream of this position
                if read1strand == '+':
                    downstreamrefpos = refpos + 1
                    try:
                        downstreamrefseq = mergedalignedpairs[downstreamrefpos][2].upper()
                        downstreamqueryseq = mergedalignedpairs[downstreamrefpos][4].upper()
                    except KeyError:
                        downstreamrefseq, downstreamqueryseq = None, None
                elif read1strand == '-':
                    downstreamrefpos = refpos - 1
                    try:
                        downstreamrefseq = mergedalignedpairs[downstreamrefpos][2].upper()
                        downstreamqueryseq = mergedalignedpairs[downstreamrefpos][4].upper()
                    except KeyError:
                        downstreamrefseq, downstreamqueryseq = None, None
                    downstreamrefseq = revcomp(downstreamrefseq)
                    downstreamqueryseq = revcomp(downstreamqueryseq)
                if downstreamrefseq == 'G':
                    conv2 = 'ng_xg'
                    if conv in ['a_x', 't_x', 'c_x']:
                        conv = None

            #Add conv(s) to dictionary
            if r2quality >= minPhred and onlyoverlap == False:
                if conv in convs:
                    convs[conv] +=1
                if conv2 == 'ng_xg' and conv != 'g_x':
                    convs[conv2] +=1

        elif r1querypos != 'NA' and r2querypos != 'NA': #this position is queried by both reads
            if read1strand == '-':
                r1refseq = revcomp(r1refseq)
                r2refseq = revcomp(r2refseq)
                r1queryseq = revcomp(r1queryseq)
                r2queryseq = revcomp(r2queryseq)
            
            if r1refseq == 'N' or r2refseq == 'N' or r1refseq == 'n' or r2refseq == 'n':
                continue
            
            r1result = r1refseq.lower() + '_' + r1queryseq.lower()
            r2result = r2refseq.lower() + '_' + r2queryseq.lower()

            #Only record if r1 and r2 agree about what is going on
            if r1result == r2result:
                conv = r1refseq.lower() + '_' + r1queryseq.lower()
                if r1queryseq == 'X':
                    #Check if there is a reference G downstream of this position
                    if read1strand == '+':
                        downstreamrefpos = refpos + 1
                        try:
                            downstreamrefseq = mergedalignedpairs[downstreamrefpos][2].upper()
                            downstreamqueryseq = mergedalignedpairs[downstreamrefpos][4].upper()
                        except KeyError:
                            downstreamrefseq, downstreamqueryseq = None, None
                    elif read1strand == '-':
                        downstreamrefpos = refpos - 1
                        try:
                            downstreamrefseq = mergedalignedpairs[downstreamrefpos][2].upper()
                            downstreamqueryseq = mergedalignedpairs[downstreamrefpos][4].upper()
                        except KeyError:
                            downstreamrefseq, downstreamqueryseq = None, None
                        downstreamrefseq = revcomp(downstreamrefseq)
                        downstreamqueryseq = revcomp(downstreamqueryseq)
                    if downstreamrefseq == 'G':
                        conv2 = 'ng_xg'
                        if conv in ['a_x', 't_x', 'c_x']:
                            conv = None

                #Add conv(s) to dictionary
                #Only do conv2 (ng_xg) if conv is not g_x
                if r1quality >= minPhred and r2quality >= minPhred:
                    if conv in convs:
                        convs[conv] +=1
                    if conv2 == 'ng_xg' and conv != 'g_x':
                        convs[conv2] +=1

        elif r1querypos == 'NA' and r2querypos == 'NA': #if we are using only read1 or read2, it's possible for this position in both reads to be NA
            continue

    #Does the number of t_c conversions meet our threshold?
    if convs['t_c'] >= nConv:
        pass
    elif convs['t_c'] < nConv:
        convs['t_c'] = 0

    return convs

def writeOutput(convs, outfile):
    #Write conv dict in text output table
    #<oligoname> <readcount> <conversion counts>
    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n']

    headerlist = ['oligo', 'readcount'] + possibleconvs
    with open(outfile, 'w') as outfh:
        outfh.write(('\t').join(headerlist) + '\n')
        for oligo in convs:
            readcount = convs[oligo][0]
            outfh.write(oligo + '\t' + str(readcount) + '\t')
            for possibleconv in possibleconvs:
                v = str(convs[oligo][1][possibleconv])
                if possibleconv != 't_n':
                    outfh.write(v + '\t')
                elif possibleconv == 't_n':
                    outfh.write(v + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Count mismatches in MPRA sequencing data.')
    parser.add_argument('--bam', type = str, help = 'Alignment file. Ideally from bowtie2.')
    parser.add_argument('--onlyConsiderOverlap', action='store_true',
                        help='Only consider conversions seen in both reads of a read pair? Only possible with paired end data.')
    parser.add_argument('--use_read1', action='store_true',
                        help='Use read1 when looking for conversions? Only useful with paired end data.')
    parser.add_argument('--use_read2', action='store_true',
                        help='Use read2 when looking for conversions? Only useful with paired end data.')
    parser.add_argument('--minMappingQual', type=int,
                        help='Minimum mapping quality for a read to be considered in conversion counting. bowtie2 unique mappers have MAPQ >=2.', required=True)
    parser.add_argument('--minPhred', type = int, help = 'Minimum phred score for a nucleotide to be considered.')
    parser.add_argument(
        '--nConv', type=int, help='Minimum number of required T->C conversions in a read pair in order for those conversions to be counted. Default is 1.', default=1)
    parser.add_argument('--output', type = str, help = 'Output file.')
    args = parser.parse_args()

    convs = iteratereads_pairedend(args.bam, args.onlyConsiderOverlap, args.use_read1, args.use_read2, args.nConv, args.minMappingQual, args.minPhred)
    writeOutput(convs, args.output)


    
