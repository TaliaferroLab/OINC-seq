#This file describes a method for trimming adapters from reads produced by libraries made using Quantseq 3' mRNA-seq V2 Library Prep Kit FWD.
#This kit is the recommended one to use with transcriptome-wide OINC-seq experiments.
#It requires that cutadapt (https://cutadapt.readthedocs.io/) be in the user's path.

#Quantseq samples
#Read1, trim AAAAAAAAAAAAAAAAAAAA from 3' end
#Read2, trim AGATCGGAAGAGCGTCGTGTAGGGAAAGACGGTA from 3' end and TTTTTTTTTTTTTTTTTTTT from 5' end

#QUANTSEQ STRATEGY
#Step1: get rid of UMI (first 6 nt of read 1) and trim 3' adapter off read 1 (-u 6 ; -a AAAAAAAAAAAAAAAAAAAA) [make temporary output files]
#Step2: Trim 5' adapter of read 2 (TTTTTTTTTTTTTTTTTTTT) [make temporary output files]
#Step3: Try to trim 3' adapter of read 2 (AGATCGGAAGAGCGTCGTGTAGGGAAAGACGGTA). Write untrimmed reads (these are done). [make temporary outfile for trimmed reads]
#Step4: For reads that did have 3' adapter on read 2, remove the last 6 bases on read 2 (UMI). These are now done too.
#Step5: Combine trimmed reads from step4 with untrimmed reads from step3.

readdir='<dir>' #Directory containing all fastq read files, eg file1_1.fq.gz, file1_2.fq.gz, file2_1.fq.gz, etc.
readfiles=( 'file1' 'file2' 'file3' ) #list of file names
outnames=( 'sample1' 'sample2' 'sample3' ) #list of final file names for trimmed reads. In this example, after trimming, file1_1.fq.gz will be named as sample1_1.fq.gz.


for i in ${!readfiles[@]}
do
	read1=${readdir}${readfiles[$i]}_1.fq.gz
	read2=${readdir}${readfiles[$i]}_2.fq.gz
	outread1s1=${readdir}${outnames[$i]}_1.temp.step1.fq.gz
	outread2s1=${readdir}${outnames[$i]}_2.temp.step1.fq.gz
	statsouts1=${readdir}${outnames[$i]}.cutadaptstats.step1.txt

	outread1s2=${readdir}${outnames[$i]}_1.temp.step2.fq.gz
    outread2s2=${readdir}${outnames[$i]}_2.temp.step2.fq.gz
    statsouts2=${readdir}${outnames[$i]}.cutadaptstats.step2.txt

	outread1s3=${readdir}${outnames[$i]}_1.temp.step3.fq.gz
    outread2s3=${readdir}${outnames[$i]}_2.temp.step3.fq.gz
    statsouts3=${readdir}${outnames[$i]}.cutadaptstats.step3.txt
	untrimmedr1s3=${readdir}${outnames[$i]}_1.untrimmed.step3.fq.gz
	untrimmedr2s3=${readdir}${outnames[$i]}_2.untrimmed.step3.fq.gz

	outread1s4=${readdir}${outnames[$i]}_1.temp.step4.fq.gz
    outread2s4=${readdir}${outnames[$i]}_2.temp.step4.fq.gz
    statsouts4=${readdir}${outnames[$i]}.cutadaptstats.step4.txt
	
	finaloutread1=${readdir}${outnames[$i]}_1.trimmed.fq.gz
	finaloutread2=${readdir}${outnames[$i]}_2.trimmed.fq.gz

	#Step1
	cutadapt -u 6 -U 0 -a AAAAAAAAAAAAAAAAAAAA --minimum-length 25 -j 8 -o ${outread1s1} -p ${outread2s1} ${read1} ${read2} > ${statsouts1}

	#Step2
	cutadapt -G TTTTTTTTTTTTTTTTTTTT --minimum-length 25 -j 8 -o ${outread1s2} -p ${outread2s2} ${outread1s1} ${outread2s1} > ${statsouts2}

	#Step3
	cutadapt -A AGATCGGAAGAGCGTCGTGTAGGGAAAGACGGTA --minimum-length 25 --untrimmed-output ${untrimmedr1s3} --untrimmed-paired-output ${untrimmedr2s3} -o ${outread1s3} -p ${outread2s3} ${outread1s2} ${outread2s2} > ${statsouts3}

	#Step4
	cutadapt -U -6 --minimum-length 25 -j 8 -o ${outread1s4} -p ${outread2s4} ${outread1s3} ${outread2s3} > ${statsouts4}

	#Step5
	cat ${untrimmedr1s3} ${outread1s4} > ${finaloutread1}
	cat ${untrimmedr2s3} ${outread2s4} > ${finaloutread2}

	rm *temp*fq.gz
	rm *untrimmed*

done