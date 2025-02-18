import os
import subprocess
import sys
import shutil
import argparse

#Given a pair of read files, align reads using STAR and quantify/align reads using salmon.
#This will make a STAR-produced bam (for pigpen mutation calling) and a salmon-produced bam (for read assignment).
#It will then run postmaster to append transcript assignments to the salmon-produced bam.

#This is going to take in gzipped fastqs, a directory containing the STAR index for this genome, and a directory containing the salmon index for this genome.

#Reads are aligned to the genome using STAR. This bam file will be used for mutation calling. 
#In this alignment, we allow multiple mapping reads, but only report the best alignment. 

#Reads are then quantified using salmon, where a separate transcriptome-oriented bam is written.
#Postmaster then takes this bam and adds posterior probabilities for transcript assignments.

#alignAndQuant.py only gives uniquely aligned reads to salmon. alignAndQuant2.py gives all reads to salmon.

#When runSTAR(), bamtofastq(), runSalmon(), and runPostmaster() are run in succession, the output is a directory called <samplename>. 
#In this directory, the STAR output is <samplename>Aligned.sortedByCoord.out.bam in STAR/,
#the salmon output is <samplename>.quant.sf and <samplename>.salmon.bam in salmon/,
#and the postmaster output is <samplename>.postmaster.bam in postmaster/

#Requires STAR, salmon(>= 1.9.0), and postmaster be in user's PATH.

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

    if reads2:
        command = ['STAR', '--runMode', 'alignReads', '--runThreadN', nthreads, '--genomeLoad', 'NoSharedMemory', '--genomeDir', STARindex, '--readFilesIn', reads1, reads2, '--readFilesCommand',
                   'zcat', '--outFileNamePrefix', prefix, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMstrandField', 'intronMotif', '--outSAMmultNmax', '1', '--outSAMattributes', 'MD', 'NH']

    elif not reads2:
        command = ['STAR', '--runMode', 'alignReads', '--runThreadN', nthreads, '--genomeLoad', 'NoSharedMemory', '--genomeDir', STARindex, '--readFilesIn', reads1, '--readFilesCommand',
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


def bamtofastq(samplename, nthreads, reads2):
    #Given a bam file of uniquely aligned reads (produced from runSTAR), rederive these reads as fastq in preparation for submission to salmon
    #This function isn't needed anymore as we will align all reads.
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
    r1file = samplename + '.unique.r1.fq.gz'
    r2file = samplename + '.unique.r2.fq.gz'
    print('Writing fastq file of uniquely aligned reads for {0}...'.format(samplename))
    if reads2:
        command = ['samtools', 'fastq', '--threads', nthreads, '-1', r1file, '-2', r2file, '-0', '/dev/null', '-s', '/dev/null', '-n', sortedbam]
    elif not reads2:
        command = ['samtools', 'fastq', '--threads', nthreads, '-0', r1file, '-n', sortedbam]
    subprocess.call(command)
    print('Done writing fastq files for {0}!'.format(samplename))

    os.remove(sortedbam)


def runSalmon(reads1, reads2, nthreads, salmonindex, samplename):

    if not os.path.exists('salmon'):
        os.mkdir('salmon')

    idx = os.path.abspath(salmonindex)
    r1 = os.path.abspath(reads1)
    if reads2:
        r2 = os.path.abspath(reads2)

    os.chdir('salmon')

    if reads2:
        command = ['salmon', 'quant', '--libType', 'A', '-p', nthreads, '--seqBias', '--gcBias',
               '--validateMappings', '-1', r1, '-2', r2, '-o', samplename, '--index', idx, '--writeMappings={0}.salmon.bam'.format(samplename), '--writeQualities']
    elif not reads2:
        command = ['salmon', 'quant', '--libType', 'A', '-p', nthreads, '--seqBias',
                   '--validateMappings', '-r', r1, '-o', samplename, '--index', idx, '--writeMappings={0}.salmon.bam'.format(samplename), '--writeQualities']

    print('Running salmon for {0}...'.format(samplename))

    subprocess.run(command)

    #Move output
    outputdir = os.path.join(os.getcwd(), samplename)
    quantfile = os.path.join(outputdir, 'quant.sf')
    movedquantfile = os.path.join(os.getcwd(), '{0}.quant.sf'.format(samplename))
    os.rename(quantfile, movedquantfile)


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
    parser.add_argument('--reversereads', type = str, help = 'Reverse reads. Gzipped fastq. Do not supply if using single end reads.')
    parser.add_argument('--nthreads', type = str, help = 'Number of threads to use for alignment and quantification.')
    parser.add_argument('--STARindex', type = str, help = 'STAR index directory.')
    parser.add_argument('--salmonindex', type = str, help = 'Salmon index directory.')
    parser.add_argument('--samplename', type = str, help = 'Sample name. Will be appended to output files.')
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

    wd = os.path.abspath(os.getcwd())
    sampledir = os.path.join(wd, samplename)
    if os.path.exists(sampledir) and os.path.isdir(sampledir):
        shutil.rmtree(sampledir)
    os.mkdir(sampledir)
    os.chdir(sampledir)

    runSTAR(r1, r2, nthreads, STARindex, samplename)
    runSalmon(r1, r2, nthreads, salmonindex, samplename)
    os.chdir(sampledir)
    runPostmaster(samplename, nthreads)

    
