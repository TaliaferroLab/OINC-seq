import os
import subprocess
import sys
import shutil
import argparse
import pysam

'''
Given a pair of read files, align reads using STAR, deduplicate reads by UMI, and quantify reads using salmon.
This will make a STAR-produced bam (for pigpen mutation calling) 
This bam will be deduplicated with UMI-tools then passed to salmon(for read assignment).
It will then run postmaster to append transcript assignments to the salmon-produced bam.

This is going to take in gzipped fastqs with UMIs extracted,
a directory containing the STAR index for this genome, and a directory containing the salmon index for this genome.

This means that, in addition to any adapter trimming, ***reads must have been first processed with umi_tools extract***.
For quantseq libraries, this corresponds to the first 6 nt of read 1.

Reads are aligned to the genome using STAR. This bam file will be used for mutation calling. 
In this alignment, we allow multiple mapping reads, but only report the best alignment. 
This bam will then be deduplicated based on UMI and alignment position.

Reads are then quantified using salmon, where a separate transcriptome-oriented bam is written.
Postmaster then takes this bam and adds posterior probabilities for transcript assignments.

When runSTAR(), runSalmon(), and runPostmaster() are run in succession, the output is a directory called <samplename>. 
In this directory, the STAR output is <samplename>.Aligned.sortedByCoord.out.bam in STAR/,
the salmon output is <samplename>.quant.sf and <samplename>.salmon.bam in salmon/,
and the postmaster output is <samplename>.postmaster.bam in postmaster/

Requires STAR, salmon(>= 1.9.0), and postmaster be in user's PATH.
'''

def runSTAR(reads1, reads2, nthreads, STARindex, samplename):
    if not os.path.exists('STAR'):
        os.mkdir('STAR')

    cwd = os.getcwd()
    outdir = os.path.join(cwd, 'STAR')

    #Clean output directory if it already exists
    if os.path.exists(outdir) and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    os.mkdir(outdir)
    prefix = outdir + '/' + samplename

    command = ['STAR', '--runMode', 'alignReads', '--runThreadN', nthreads, '--genomeLoad', 'NoSharedMemory', '--genomeDir', STARindex, '--readFilesIn', reads1, reads2, '--readFilesCommand',
               'zcat', '--outFileNamePrefix', prefix, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMstrandField', 'intronMotif', '--outSAMmultNmax', '1', '--outSAMattributes', 'MD', 'NH']

    print('Running STAR for {0}...'.format(samplename))

    subprocess.run(command)

    #make index
    bam = os.path.join(outdir, samplename + 'Aligned.sortedByCoord.out.bam')
    bamindex = bam + '.bai'
    if not os.path.exists(bamindex):
        indexCMD = 'samtools index ' + bam
        index = subprocess.Popen(indexCMD, shell=True)
        index.wait()

    print('Finished STAR for {0}!'.format(samplename))


def filterbam(samplename, maxmap):
    #Take a bam and filter it, only keeping reads that map to <= maxmap locations using NH:i tag
    #For some reason whether STAR uses --outFilterMultiMapNmax is unpredictably variable, so we will do it this way.
    cwd = os.getcwd()
    outdir = os.path.join(cwd, 'STAR')
    maxmap = int(maxmap)
    inbam = os.path.join(outdir, samplename + 'Aligned.sortedByCoord.out.bam')
    outbam = os.path.join(outdir, samplename +
                          'Aligned.sortedByCoord.multifiltered.out.bam')

    print('Removing reads with > {0} alignments...'.format(maxmap))

    with pysam.AlignmentFile(inbam, 'rb') as infh, pysam.AlignmentFile(outbam, 'wb', template=infh) as outfh:
        readcount = 0
        filteredreadcount = 0
        for read in infh.fetch(until_eof=True):
            readcount += 1
            nh = read.get_tag('NH')
            if nh <= maxmap:
                filteredreadcount += 1
                outfh.write(read)

    #Remove unfiltered bam and its index
    os.remove(inbam)
    os.remove(inbam + '.bai')
    #Rename filtered bam so that it has the same name as the original
    #This helps later when pipgen is looking for bams with certain expected names
    os.rename(outbam, inbam)
    #index filtered bam
    bamindex = inbam + '.bai'
    indexCMD = 'samtools index ' + inbam
    index = subprocess.Popen(indexCMD, shell=True)
    index.wait()

    filteredpct = round((filteredreadcount / readcount) * 100, 3)

    print('Looked through {0} reads. {1} ({2}%) had {3} or fewer alignments.'.format(
        readcount, filteredreadcount, filteredpct, maxmap))


def runDedup(samplename, nthreads):
    STARbam = os.path.join(os.getcwd(), 'STAR', '{0}Aligned.sortedByCoord.out.bam'.format(samplename))
    dedupbam = os.path.join(os.getcwd(), 'STAR', '{0}.dedup.bam'.format(samplename))
    if args.libType == "LEXO":
        command = ['umi_tools', 'dedup', '-I', STARbam, '--paired', '-S', dedupbam]
    elif args.libType == "SA":
        command = ['umi_tools', 'dedup', '-I', STARbam, '--paired', '--method=unique', '-S', dedupbam]
    else:
        print('LibType must be either "LEXO" or "SA".')
    print('Running deduplication for {0}...'.format(samplename))

    subprocess.run(command)

    command = ['samtools', 'index', dedupbam]
    subprocess.run(command)

    #We don't need the STAR alignment file anymore, and it's pretty big
    #Rename to the old name so downstream code finds the bams it's looking for
    os.rename(dedupbam, STARbam)
    #Reindex deduplicated bam
    bamindex = STARbam + '.bai'
    indexCMD = 'samtools index ' + STARbam
    index = subprocess.Popen(indexCMD, shell=True)
    index.wait()

    print('Finished deduplicating {0}!'.format(samplename))


def bamtofastq(samplename, nthreads, dedup, reads2):
    #Given a bam file of uniquely aligned reads (produced from runSTAR), rederive these reads as fastq in preparation for submission to salmon
    if not os.path.exists('STAR'):
        os.mkdir('STAR')

    cwd = os.getcwd()
    outdir = os.path.join(cwd, 'STAR')
    inbam = os.path.join(outdir, samplename + 'Aligned.sortedByCoord.out.bam')
    sortedbam = os.path.join(outdir, 'temp.namesort.bam')

    #First sort bam file by readname
    print('Sorting bam file by read name...')
    command = ['samtools', 'collate', '--threads', nthreads, '-u', '-o', sortedbam, inbam]
    subprocess.call(command)
    print('Done!')

    #Now derive fastq
    r1file = samplename + '.aligned.r1.fq.gz'
    r2file = samplename + '.aligned.r2.fq.gz'
    if dedup:
        print('Writing fastq file of deduplicated reads for {0}...'.format(samplename))
    elif not dedup:
        print('Writing fastq file of aligned reads for {0}...'.format(samplename))
    if reads2:
        command = ['samtools', 'fastq', '--threads', nthreads, '-1', r1file, '-2', r2file, '-0', '/dev/null', '-s', '/dev/null', '-n', sortedbam]
    elif not reads2:
        command = ['samtools', 'fastq', '--threads', nthreads, '-0', r1file, '-n', sortedbam]
    subprocess.call(command)
    print('Done writing fastq files for {0}!'.format(samplename))
    
    os.remove(sortedbam)


def runSalmon(reads1, reads2, nthreads, salmonindex, samplename):
    #Take in those deduplicated reads and quantify transcript abundance with them using salmon.

    if not os.path.exists('salmon'):
        os.mkdir('salmon')

    idx = os.path.abspath(salmonindex)
    r1 = os.path.abspath(reads1)
    r2 = os.path.abspath(reads2)

    os.chdir('salmon')

    command = ['salmon', 'quant', '--libType', 'A', '-p', nthreads, '--seqBias', '--gcBias',
               '--validateMappings', '-1', r1, '-2', r2, '-o', samplename, '--index', idx, '--writeMappings={0}.salmon.bam'.format(samplename), '--writeQualities']

    print('Running salmon for {0}...'.format(samplename))

    subprocess.call(command)

    #Move output
    outputdir = os.path.join(os.getcwd(), samplename)
    quantfile = os.path.join(outputdir, 'quant.sf')
    movedquantfile = os.path.join(os.getcwd(), '{0}.quant.sf'.format(samplename))
    os.rename(quantfile, movedquantfile)

    #Remove uniquely aligning read files
    os.remove(r1)
    if reads2:
        os.remove(r2)

    print('Finished salmon for {0}!'.format(samplename))


def runPostmaster(samplename, nthreads):
    if not os.path.exists('postmaster'):
        os.mkdir('postmaster')

    salmonquant = os.path.join(os.getcwd(), 'salmon', '{0}.quant.sf'.format(samplename))
    salmonbam = os.path.join(os.getcwd(), 'salmon', '{0}.salmon.bam'.format(samplename))

    os.chdir('postmaster')
    outputfile = os.path.join(os.getcwd(), '{0}.postmaster.bam'.format(samplename))

    print('Running postmaster for {0}...'.format(samplename))
    command = ['postmaster', '--num-threads', nthreads, '--quant', salmonquant, '--alignments', salmonbam, '--output', outputfile]
    subprocess.call(command)

    #Sort and index bam
    with open(outputfile + '.sort', 'w') as sortedfh:
        command = ['samtools', 'sort', '-@', nthreads, outputfile]
        subprocess.run(command, stdout = sortedfh)
    os.rename(outputfile + '.sort', outputfile)

    command = ['samtools', 'index', outputfile]
    subprocess.run(command)

    #We don't need the salmon alignment file anymore, and it's pretty big
    os.remove(salmonbam)

    print('Finished postmaster for {0}!'.format(samplename))


def addMD(samplename, reffasta, nthreads):
    inputbam = os.path.join(os.getcwd(), 'postmaster', '{0}.postmaster.bam'.format(samplename))
    command = ['samtools', 'calmd', '-b', '--threads', nthreads, inputbam, reffasta]

    print('Adding MD tags to {0}.postmaster.md.bam...'.format(samplename))
    with open(samplename + '.postmaster.md.bam', 'w') as outfile:
        subprocess.run(command, stdout = outfile)
    print('Finished adding MD tags to {0}.postmaster.md.bam!'.format(samplename))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Align and quantify reads using STAR, salmon, and postmaster in preparation for analysis with PIGPEN.')
    parser.add_argument('--forwardreads', type = str, help = 'Forward reads. Gzipped fastq.')
    parser.add_argument('--reversereads', type = str, help = 'Reverse reads. Gzipped fastq.')
    parser.add_argument('--nthreads', type = str, help = 'Number of threads to use for alignment and quantification.')
    parser.add_argument('--STARindex', type = str, help = 'STAR index directory.')
    parser.add_argument('--salmonindex', type = str, help = 'Salmon index directory.')
    parser.add_argument('--samplename', type = str, help = 'Sample name. Will be appended to output files.')
    parser.add_argument('--dedupUMI', action = 'store_true', help = 'Deduplicate UMIs? requires UMI extract.')
    parser.add_argument('--libType', type = str, help = 'Library Type, either "LEXO" or "SA"')
    parser.add_argument(
        '--maxmap', type=int, help='Maximum number of allowable alignments for a read.')
    args = parser.parse_args()

    r1 = os.path.abspath(args.forwardreads)
    if args.reversereads:
        r2 = os.path.abspath(args.reversereads)
    elif not args.reversereads:
        r2 = None
    STARindex = os.path.abspath(args.STARindex)
    salmonindex = os.path.abspath(args.salmonindex)
    samplename = args.samplename
    nthreads = args.nthreads
    maxmap = args.maxmap

    wd = os.path.abspath(os.getcwd())
    sampledir = os.path.join(wd, samplename)
    if os.path.exists(sampledir) and os.path.isdir(sampledir):
        shutil.rmtree(sampledir)
    os.mkdir(sampledir)
    os.chdir(sampledir)


    runSTAR(r1, r2, nthreads, STARindex, samplename)
    filterbam(samplename, maxmap)
    if args.dedupUMI:
        runDedup(samplename, nthreads)

    #aligned read files
    alignedr1 = samplename + '.aligned.r1.fq.gz'
    if args.reversereads:
        alignedr2 = samplename + '.aligned.r2.fq.gz'
    elif not args.reversereads:
        alignedr2 = None
    bamtofastq(samplename, nthreads, args.dedupUMI, r2)
    runSalmon(alignedr1, alignedr2, nthreads, salmonindex, samplename)
    #Remove aligned fastqs
    os.chdir(sampledir)
    os.chdir(sampledir)
    runPostmaster(samplename, nthreads)
