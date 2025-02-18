#python >=3.6
import random
import sys
import gzip

#usage: python simulateOINCreads.py <readlength> <readdepth>

#This is the wildtype sequence of the amplicon. Paired end reads will read in from both ends.
seq = 'ACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACACCCACTCCTCCACCTTTGACGCTGGGGCTGGCATTGCCCTCAACGACCACTTTGTCAAGCTC'
#first 200 nt of above seq so that we can incorporate read overlap
seq = 'ACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGT'


#Intrinsic (i.e. cell- or RT-derived) mutation rates
mutfreqs = {'A' : {'C' : 1e-5, 'T' : 0, 'G' : 2e-4},
'C' : {'G' : 5e-5, 'T' : 1.5e-3, 'A' : 4e-4},
'G' : {'C' : 3e-5, 'T' : 2e-6, 'A' : 8e-4},
'T' : {'A' : 3e-6, 'C' : 7e-5, 'G' : 8e-6}}

#Sequencing error rate
seqerrorrate = 0.001


def revcomp(nt):
    revcompdict = {
        'G': 'C',
        'C': 'G',
        'A': 'T',
        'T': 'A',
        'N': 'N',
        'g': 'c',
        'c': 'g',
        'a': 't',
        't': 'a',
        'n': 'n'
    }

    nt_rc = revcompdict[nt]

    return nt_rc

def makecDNAseq(wtseq, mutfreqs):
    #Make the sequence of the cDNA for this read pair.
    #This is intended to be the entire fragment. We will break it up
    #into the read pairs (and incorporate sequencing errors) later.
    #build sequence nt by nt
    outseq = ''
    for nt in wtseq:
        possiblents = list(mutfreqs[nt].keys())
        possiblentfreqs = list(mutfreqs[nt].values())
        wtfreq = 1 - sum(possiblentfreqs)
        #add chance nt will be wt
        possiblents.append(nt)
        possiblentfreqs.append(wtfreq)
        
        outnt = random.choices(
            population = possiblents,
            weights = possiblentfreqs,
            k = 1
        )
        outnt = outnt[0]
        outseq += outnt

    return outseq

def addseqerrors(readseq, seqerrorrate):
    #Given a read sequence, add simulated sequencing erros
    outseq = ''
    for nt in readseq:
        #is there a sequencing error at this position?
        possiblents = list(mutfreqs[nt].keys())
        possiblentfreqs = [seqerrorrate / 3] * 3
        #add chance of not having a sequencing error
        possiblents.append(nt)
        possiblentfreqs.append(1 - seqerrorrate)
        
        outnt = random.choices(
            population = possiblents,
            weights = possiblentfreqs,
            k = 1
        )
        outnt = outnt[0]
        outseq += outnt

    return outseq

def makefastq(seq, mutfreqs, seqerrorrate, readlength, depth, outfile):
    #Given a wildtype amplicon (seq), paired end read length, desired depth, make simulated reads
    #with desired mutation and sequencing error rates
    #all quality scores are J (41)
    readlength = int(readlength)
    depth = int(depth)

    with gzip.open(outfile + '_1.fq.gz', 'wt') as read1outfh, gzip.open(outfile + '_2.fq.gz', 'wt') as read2outfh:
        readcounter = 0
        for i in range(depth):
            readcounter +=1
            if readcounter % 100000 == 0:
                print('Creating read {0}...'.format(readcounter))

            #Make fragment for this readpair
            fragseq = makecDNAseq(seq, mutfreqs)
            #Make reads from this fragment
            read1seq = fragseq[:readlength]
            read2seq = fragseq[readlength * -1:]
            #reverse complement read2
            read2seq_rc = ''
            for nt in read2seq:
                nt_rc = revcomp(nt)
                read2seq_rc += nt_rc
            read2seq_rc = read2seq_rc[::-1]

            #Add sequencing errors
            read1seq = addseqerrors(read1seq, seqerrorrate)
            read2seq = addseqerrors(read2seq_rc, seqerrorrate)

            qualityscores = 'J' * len(read1seq)
            readtitle = '@simread_' + str(readcounter)

            read1outfh.write(readtitle + '\n' + read1seq + '\n' + '+' + '\n' + qualityscores + '\n')
            read2outfh.write(readtitle + '\n' + read2seq + '\n' + '+' + '\n' + qualityscores + '\n')

    with open('simulationparams.txt', 'w') as outfh:
        outfh.write(('\t').join(['ref', 'mut', 'freq']) + '\n')
        for ref in mutfreqs:
            for mut in mutfreqs[ref]:
                freq = str(mutfreqs[ref][mut])
                outfh.write(('\t').join([ref, mut, freq]) + '\n')


makefastq(seq, mutfreqs, seqerrorrate, sys.argv[1], sys.argv[2], 'oincsimulation')



