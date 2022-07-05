import os
import subprocess
import sys
import shutil
import argparse

#Given a pair of read files, align reads using STAR and quantify/align reads using salmon.
#This will make a STAR-produced bam (for pigpen mutation calling) and a salmon-produced bam (for read assignment).
#It will then run postmaster to append transcript assignments to the salmon-produced bam.

#This is going to take in gzipped fastqs, a directory containing the STAR index for this genome, and a directory containing the salmon index for this genome.

#Reads are aligned to the genome using STAR. This bam file will be used for mutation calling. Uniquely aligning reads from this alignment are then
#written to temporary fastq files (<samplename>.unique.r1..fq.gz), which are then used for salmon and postmaster.

#When runSTAR(), bamtofastq(), runSalmon(), and runPostmaster() are run in succession, the output is a file called <samplename>.postmaster.bam in the postmaster/ 
#and <samplename>Aligned.sortedByCoord.out.bam in STAR/<samplename> and <samplename>.quant.sf in salmon/

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

    command = ['STAR', '--runMode', 'alignReads', '--runThreadN', nthreads, '--genomeLoad', 'NoSharedMemory', '--genomeDir', STARindex, '--readFilesIn', reads1, reads2, '--readFilesCommand',
               'zcat', '--outFileNamePrefix', prefix, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMstrandField', 'intronMotif', '-–outFilterMultimapNmax', '1', '--outSAMattributes', 'MD', 'NH']

    print('Running STAR for {0}...'.format(samplename))
    
    subprocess.call(command)

    print('Finished STAR for {0}!'.format(samplename))


def bamtofastq(samplename, nthreads):
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
    r1file = samplename + '.unique.r1.fq.gz'
    r2file = samplename + '.unique.r2.fq.gz'
    print('Writing fastq file of uniquely aligned reads for {0}...'.format(samplename))
    command = ['samtools', 'fastq', '--threads', nthreads, '-1', r1file, '-2', r2file, '-0', '/dev/null', '-s', '/dev/null', '-n', sortedbam]
    subprocess.call(command)
    print('Done writing fastq files for {0}!'.format(samplename))

    os.remove(sortedbam)


def runSalmon(reads1, reads2, nthreads, salmonindex, samplename):
    #Take in those uniquely aligning reads and quantify transcript abundance with them using salmon.

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
    subprocess.call(command)
    
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
    args = parser.parse_args()

    r1 = os.path.abspath(args.forwardreads)
    r2 = os.path.abspath(args.reversereads)
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

    #uniquely aligning read files
    uniquer1 = samplename + '.unique.r1.fq.gz'
    uniquer2 = samplename + '.unique.r2.fq.gz'

    runSTAR(r1, r2, nthreads, STARindex, samplename)
    bamtofastq(samplename, nthreads)
    runSalmon(uniquer1, uniquer2, nthreads, salmonindex, samplename)
    os.chdir(sampledir)
    runPostmaster(samplename, nthreads)

    
