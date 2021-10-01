#Given a bam, first filter it for reads that lie in regions we care about
#For example, if this is a 3' end sequencing library, filter it for reads
#that lie in last exons (as defined by a bed).
#The idea here is that this will (1) potentially remove a lot of iffy alignments that could
#have multiple mismatches and (2) speed up the counting of mismatches because we won't have
#to find them in every read.

#Both bam and bed must be sorted by coordinate

import pybedtools
import os
import sys
import subprocess

def intersectreads(bam, bed, chrsort):
    #Get a list of reads that overlap bed intervals
    #chrsort is a tab-delimited file of chromosomes in the order the appear in the bed/bam and their sizes
    #can be made with cut -f 1,2 genome.fa.fai
    bam = pybedtools.BedTool(bam)
    bed = pybedtools.BedTool(bed)
    readstokeep = []

    intersectingreads = bam.intersect(bed, u = True, bed = True, sorted = True, g = chrsort)
   
    for read in intersectingreads:
        read = str(read)
        read = read.strip().split('\t')
        readid = read[3]

        #bedtools gives you each read of a mate pair individually
        #samtools (later) is going to filter on the pair ID (i.e. without /1 or /2)
        #so if either read of a pair overlaps, it's going to end up in our filtered bam
        readbaseid = readid[:-2]
        readstokeep.append(readbaseid)

    #remove duplicates
    readstokeep = list(set(readstokeep))

    with open('reads.tmp', 'w') as outfh:
        for read in readstokeep:
            outfh.write(read + '\n')


def filterbam(bam):
    #First see how many alignments are in this file
    countCMD = 'samtools view -c ' + bam
    countCMD = countCMD.split(' ')
    alncount = subprocess.check_output(countCMD).decode('utf-8').strip()
    print('Filtering {0} alignments in {1}...'.format(alncount, os.path.basename(bam)))

    filteredbamBasename = os.path.basename(bam).replace('.bam', '.filtered.bam')
    filteredbamPath = os.path.join(os.path.dirname(os.path.abspath(bam)), filteredbamBasename)
    filterCMD = 'samtools view -N reads.tmp -o ' + filteredbamPath + ' ' + bam
    filt = subprocess.Popen(filterCMD, shell = True)
    filt.wait()

    countCMD = 'samtools view -c ' + filteredbamPath
    countCMD = countCMD.split(' ')
    filteredalncount = subprocess.check_output(countCMD).decode('utf-8').strip()
    filterpct = round(int(filteredalncount) / int(alncount), 3) * 100

    print('{0} of {1} alignments are in filtered bam file ({2}%)'.format(filteredalncount, alncount, filterpct))

    os.remove('reads.tmp')


if __name__ == '__main__':
    intersectreads(sys.argv[1], sys.argv[2], sys.argv[3])
    filterbam(sys.argv[1])



