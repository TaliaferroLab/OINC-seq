#We want to mask snps as they could be contributing "errors" that are not dependent on 8oxoG.
#Take bams (probably bams from control, minusDBF samples) and detect variants.
#Then mask these positions in downstream analyses.

#If there are multiple bams to detect variants from, merge their vcf files.

import subprocess
import os
import sys

#This will take in a list of bams and identify variants, creating vcf files for each
def getSNPs(bams, genomefasta, minCoverage = 20, minVarFreq = 0.02):
    if not minCoverage:
        minCoverage = 20
    if not minVarFreq:
        minVarFreq = 0.02

    #if we already made a vcf, don't make another one
    if os.path.exists('merged.vcf'):
        print('A merged vcf files already exists! Not making another one...')
        return None

    vcfFileNames = []

    for bam in bams:
        vcffn = os.path.basename(bam).replace('.bam', '.vcf')
        vcfFileNames.append(vcffn + '.gz')
        vcflogfn = vcffn + '.log'
        with open(vcffn, 'w') as vcffh, open(vcflogfn, 'w') as logfh:
            print('Piling up reads in {0}...'.format(bam))
            mpileupCMD = 'samtools mpileup -B -A -f ' + genomefasta + ' ' + bam
            mpileup = subprocess.Popen(mpileupCMD, shell = True, stdout = subprocess.PIPE, stderr = logfh)

            varscanCMD = 'varscan mpileup2snp --strand-filter 0 --output-vcf --min-BQ 20 --min-var-freq ' + str(minVarFreq) + ' --min-coverage ' + str(minCoverage) + ' --variants 1'
            
            print('Finding variants in {0}...'.format(bam))
            varscan = subprocess.Popen(varscanCMD, shell = True, stdin = mpileup.stdout, stdout = vcffh, stderr = logfh)
            varscan.wait()

            #compress
            bgzipCMD = 'bgzip -f ' + vcffn
            compress = subprocess.Popen(bgzipCMD, shell = True)
            compress.wait()


            #make index
            indexCMD = 'bcftools index -f ' + vcffn + '.gz'
            icmd = subprocess.Popen(indexCMD, shell = True)
            icmd.wait()

    #if there is more than 1 vcf file, merge them
    if len(vcfFileNames) == 1:
        with open('vcfconcat.log', 'w') as logfh:
            vcfFileNames = vcfFileNames * 2
            vcfFiles = ' '.join(vcfFileNames)
            print(vcfFiles)
            concatCMD = 'bcftools merge --force-samples -m snps -O z --output merged.vcf ' + vcfFiles
            concat = subprocess.Popen(concatCMD, shell = True, stderr = logfh)
            concat.wait()
        filetorecord = vcfFileNames[0]
    elif len(vcfFileNames) > 1:
        with open('vcfconcat.log', 'w') as logfh:
            vcfFiles = ' '.join(vcfFileNames)
            concatCMD = 'bcftools merge --force-samples -m snps -O z --output merged.vcf ' + vcfFiles
            concat = subprocess.Popen(concatCMD, shell = True, stderr = logfh)
            concat.wait()

    return vcfFileNames

#Go through each vcf, and record locations (chr:position) of SNPs so that they can be masked
#vcf coordinates are 1-based
#turn them into 0-based coordinates so that they can easily be intersected with pysam get_aligned_pairs(), which returns 0-based coordinates
def recordSNPs(mergedvcf):
    snps = {} #locations of snps {chr : set(positions)}
    with open(mergedvcf, 'r') as infh:
        for line in infh:
            if '#' in line:
                continue
            line = line.strip().split('\t')
            chrom = line[0]
            pos = int(line[1])
            if chrom not in snps:
                snps[chrom] = [pos - 1]
            elif chrom in snps:
                snps[chrom].append(pos - 1)

    #Turn into sets for faster membership searching
    for chrom in snps:
        snps[chrom] = set(snps[chrom])

    #Tell us how many snps we found
    totalsnps = 0
    for chrm in snps:
        x = len(snps[chrm])
        totalsnps += x

    print('Masking {0} variant positions.'.format(totalsnps))

    return snps


if __name__ == '__main__':
    #2 bams here just for testing
    bams = sys.argv[1:5]
    genomefasta = sys.argv[5]
    minCoverage = int(sys.argv[6])
    minVarFreq = float(sys.argv[7])

    vcfFileNames = getSNPs(bams, genomefasta, minCoverage, minVarFreq)
    recordSNPs(vcfFileNames)
