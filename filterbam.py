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
import multiprocessing as mp
from getmismatches import split_bam

def split_bed(bed):
    #Split a bed file into many bed files, one per chromosome

    #make dictionary of bed entries
    beddict = {} #{chr : [[bed line]]}
    with open(bed, 'r') as infh:
        for line in infh:
            line = line.strip().split('\t')
            chrm = line[0]
            if chrm in beddict:
                beddict[chrm].append(line)
            elif chrm not in beddict:
                beddict[chrm] = [line]

    chrms = list(beddict.keys())
    splitbeds = []

    for chrm in chrms:
        outfilebasename = chrm + '.temp.bed'
        outfilename = os.path.join(os.path.dirname(os.path.abspath(bed)), outfilebasename)
        splitbeds.append(outfilename)
        with open(outfilename, 'w') as outfh:
            bedlines = beddict[chrm]
            for line in bedlines:
                outfh.write(('\t').join(line) + '\n')

    #add path stem to splitbeds
    splitbeds = [os.path.join(os.path.dirname(os.path.abspath(bed)), f) for f in splitbeds]
    
    return splitbeds

    

    
def intersectreads(bam, bed, chrsort, writereads = False):
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

    if writereads:
        with open('reads.tmp', 'w') as outfh:
            for read in readstokeep:
                outfh.write(read + '\n')

    else:
        return readstokeep


def filterbam(bam, nproc):
    #First see how many alignments are in this file
    countCMD = 'samtools view -c ' + bam
    countCMD = countCMD.split(' ')
    alncount = subprocess.check_output(countCMD).decode('utf-8').strip()
    print('Filtering {0} alignments in {1}...'.format(alncount, os.path.basename(bam)))

    filteredbamBasename = os.path.basename(bam).replace('.bam', '.filtered.bam')
    filteredbamPath = os.path.join(os.path.dirname(os.path.abspath(bam)), filteredbamBasename)
    filterCMD = 'samtools view -@ ' + str(nproc) + ' -N reads.tmp -o ' + filteredbamPath + ' ' + bam
    filt = subprocess.Popen(filterCMD, shell = True)
    filt.wait()

    countCMD = 'samtools view -c ' + filteredbamPath
    countCMD = countCMD.split(' ')
    filteredalncount = subprocess.check_output(countCMD).decode('utf-8').strip()
    filterpct = round(int(filteredalncount) / int(alncount), 3) * 100

    print('{0} of {1} alignments are in filtered bam file ({2}%)'.format(filteredalncount, alncount, filterpct))

    os.remove('reads.tmp')

    return filteredbamPath

def intersectreads_multiprocess(bam, bed, chrsort, nproc):
    #intersect reads (calling intersectreads() above), but with bams and beds split by chromosome
    #if there's only one processor, easier to use intersectreads() directly
    splitbeds = split_bed(bed) #splits bed and returns filenames of split beds
    splitbams = split_bam(bam, int(nproc)) #splits bam and return filenames of split bams
    
    #index bams
    for f in splitbams:
        indexCMD = 'samtools index ' + f
        index = subprocess.Popen(indexCMD, shell = True)
        index.wait()

    pool = mp.Pool(processes = int(nproc))
    argslist = []
    for splitbed in splitbeds:
        chrm = os.path.basename(splitbed).replace('.temp.bed', '')
        bamid = chrm + '_SPLIT.bam'
        splitbam = [x for x in splitbams if bamid in x]
        if len(splitbam) != 1:
            print('ERROR: can\'t find split bam associated with chromosome {0}!.'.format(chrm))
            sys.exit()
        else:
            splitbam = splitbam[0]

        argslist.append((splitbam, splitbed, chrsort, False))

    print('Finding reads that intersect with bed intervals using {0} processors.'.format(nproc))
    results = pool.starmap(intersectreads, argslist)
    #write reads to file
    with open('reads.tmp', 'w') as outfh:
        for x in results:
            for read in x:
                outfh.write(read + '\n')

    #clean up
    for f in splitbeds:
        os.remove(f)
    for f in splitbams:
        os.remove(f)
        os.remove(f + '.bai')

    







if __name__ == '__main__':
    intersectreads_multiprocess(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    #intersectreads(sys.argv[1], sys.argv[2], sys.argv[3])
    filterbam(sys.argv[1], sys.argv[4])

    #split_bam_samtools(sys.argv[1], sys.argv[2])



