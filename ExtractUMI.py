import os
import subprocess
import sys
import shutil
import argparse
'''
Given a pair of read files, extract UMIs matching a given pattern
This step is required for any downstream UMI deduplication
Requires UMI_tools 
Usage:
python -u ~/Projects/OINC_seq/test/py2test/ExtractUMI.py --forwardreads ./ATP5MC1_Rep1_pDBF.10M.R1.fq.gz,./ATP5MC1_Rep1_mDBF.10M.R1.fq.gz --reversereads ./ATP5MC1_Rep1_pDBF.10M.R2.fq.gz,./ATP5MC1_Rep1_mDBF.10M.R2.fq.gz --samplename pDBF10M,mDBF10M --lib_type LEXO

'''
def runExtract(r1, r2, samplename, lib_type):
    if not os.path.exists('UMI_fastq'):
        os.mkdir('UMI_fastq')

    cwd = os.getcwd()
    outdir = os.path.join(cwd, 'UMI_fastq')
   
    r1 = r1.split(",")
    r2 = r2.split(",")
    samplename = samplename.split(",")
    
    for idx, sample in enumerate(samplename):
        
        reads1 = r1[idx]
        reads2 = r2[idx]
        output1 = outdir + '/' + sample + '.R1.fq.gz'
        output2 = outdir + '/' + sample + '.R2.fq.gz'
    
        if lib_type == "LEXO":
            command = ["umi_tools", "extract", "-I", reads1, '--bc-pattern=NNNNNN', '--read2-in={0}'.format(reads2), '--stdout={0}'.format(output1),'--read2-out={0}'.format(output2)]
        elif lib_type == "SA":
            command = ["umi_tools", "extract", "-I", reads2, '--bc-pattern=NNNNNNNNNNNN', '--read2-in={0}'.format(reads2), '--stdout={0}'.format(output1),'--read2-out={0}'.format(output2)]
        else:
            print('--lib_type must be either "LEXO" or "SA"')
            sys.exit()
    
        print('Extracting UMIs for {0}...'.format(sample))
        subprocess.run(command)
        print('Finished Extracting UMIs for {0}!'.format(sample))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Extract UMIs using umi-tools in preparation for analysis with AlignUMIquant.')
    parser.add_argument('--forwardreads', type = str, help = 'Forward reads. Gzipped fastq.', required = True)
    parser.add_argument('--reversereads', type = str, help = 'Reverse reads. Gzipped fastq.', required = True)
    parser.add_argument('--samplename', type = str, help = 'Sample name. Will be appended to output files.', required = True)
    parser.add_argument('--lib_type', type = str, help = 'Library type. Either "LEXO" or "SA"', required = True)
    args = parser.parse_args()
    
    r1 = os.path.abspath(args.forwardreads)
    r2 = os.path.abspath(args.reversereads)    
    samplename = args.samplename
    lib_type = args.lib_type
    
    if args.lib_type not in ["LEXO", "SA"]:
        print('--lib_type must be either "LEXO" or "SA"')
        sys.exit()

    runExtract(r1, r2, samplename, lib_type)

