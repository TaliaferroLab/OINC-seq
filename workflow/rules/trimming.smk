"""
Rules for trimming NGS reads with cutadapt
(http://cutadapt.readthedocs.org/en/latest/guide.html#illumina-truseq)
For usage, include this in your workflow.

Quantseq samples
Read1, trim AAAAAAAAAAAAAAAAAAAA from 3' end
Read2, trim AGATCGGAAGAGCGTCGTGTAGGGAAAGACGGTA from 3' end and TTTTTTTTTTTTTTTTTTTT from 5' end

NEW QUANTSEQ STRATEGY
Step1: get rid of random hex bindng site (first 6 nt of read 1) and trim 3' adapter off read 1 (-u 6 ; -a AAAAAAAAAAAAAAAAAAAA) [make temporary output files]
Step2: Trim 5' adapter of read 2 (TTTTTTTTTTTTTTTTTTTT) [make temporary output files]
Step3: Try to trim 3' adapter of read 2 (AGATCGGAAGAGCGTCGTGTAGGGAAAGACGGTA). Write untrimmed reads (these are done). [make temporary outfile for trimmed reads]
Step4: For reads that did have 3' adapter on read 2, remove the last 6 bases on read 2 (UMI). These are now done too.
Step5: Combine trimmed reads from step4 with untrimmed reads from step3.

SINGLE AMPLICON (SA) STRATEGY
Step1: get rid of primer binding sites (20nt) + internal barcode (6nt)
"""
if not config["libtype"]:
    print("libtype must be included in the config file")

elif config["libtype"]=="LEXO":
    rule cutadapt_LEXO:
        """Trims given paired-end reads with given parameters"""
        input:
            read1=expand("UMI_fastq/{sample}_R1.fq.gz", sample=config["samples"]),
            read2=expand("UMI_fastq/{sample}_R2.fq.gz", sample=config["samples"]),        
        output:
            outread1s1=temp(expand("cutadapt/{sample}_1.temp.step1.fq.gz", sample=config["samples"])),
            outread2s1=temp(expand("cutadapt/{sample}_2.temp.step1.fq.gz", sample=config["samples"])),
            statsouts1=expand("cutadapt/{sample}.cutadaptstats.step1.txt", sample=config["samples"]),
            outread1s2=temp(expand("cutadapt/{sample}_1.temp.step2.fq.gz", sample=config["samples"])),
            outread2s2=temp(expand("cutadapt/{sample}_2.temp.step2.fq.gz", sample=config["samples"])),
            statsouts2=expand("cutadapt/{sample}.cutadaptstats.step2.txt", sample=config["samples"]),
            outread1s3=temp(expand("cutadapt/{sample}_1.temp.step3.fq.gz", sample=config["samples"])),
            outread2s3=temp(expand("cutadapt/{sample}_2.temp.step3.fq.gz", sample=config["samples"])),
            statsouts3=expand("cutadapt/{sample}.cutadaptstats.step3.txt", sample=config["samples"]),
            untrimmedr1s3=temp(expand("cutadapt/{sample}_1.untrimmed.step3.fq.gz", sample=config["samples"])),
            untrimmedr2s3=temp(expand("cutadapt/{sample}_2.untrimmed.step3.fq.gz", sample=config["samples"])),
            outread1s4=temp(expand("cutadapt/{sample}_1.temp.step4.fq.gz", sample=config["samples"])),
            outread2s4=temp(expand("cutadapt/{sample}_2.temp.step4.fq.gz", sample=config["samples"])),
            statsouts4=expand("cutadapt/{sample}.cutadaptstats.step4.txt", sample=config["samples"]), 
            finaloutread1=expand("cutadapt/{sample}_1.trimmed.fq.gz", sample=config["samples"]),
            finaloutread2=expand("cutadapt/{sample}_2.trimmed.fq.gz", sample=config["samples"]),               
        threads: config["threads"]
        run:
            for Ir1,Ir2,O1r1,O1r2,S1,O2r1,O2r2,S2,O3r1,O3r2,S3,U3r1,U3r2,O4r1,O4r2,S4,OFr1,OFr2 in zip(input.read1,input.read2,output.outread1s1,output.outread2s1,output.statsouts1,output.outread1s2,output.outread2s2,output.statsouts2,output.outread1s3,output.outread2s3,output.statsouts3,output.untrimmedr1s3,output.untrimmedr2s3,output.outread1s4,output.outread2s4,output.statsouts4,output.finaloutread1,output.finaloutread2):
                #Step1
                shell(
                    "cutadapt -u 6 -U 0 -a AAAAAAAAAAAAAAAAAAAA --minimum-length 25 "
                    "-j {threads} -o {O1r1} -p {O1r2} {Ir1} {Ir2} > {S1}"
                )
        
                #Step2
                shell(
                    "cutadapt -G TTTTTTTTTTTTTTTTTTTT --minimum-length 25 "
                    "-j {threads} -o {O2r1} -p {O2r2} {O1r1} {O1r2} > {S2}"
                )
        
                #Step3
                shell(
                    "cutadapt -A AGATCGGAAGAGCGTCGTGTAGGGAAAGACGGTA --minimum-length 25 "
                    "-j {threads} --untrimmed-output {U3r1} "
                    "--untrimmed-paired-output {U3r2} -o {O3r1} -p {O3r2} "
                    "{O2r1} {O2r2} > {S3}"
                )
        
                #Step4
                shell(
                    "cutadapt -U -6 --minimum-length 25 -j {threads} -o {O4r1} "
                    "-p {O4r2} {O3r1} {O3r2} > {S4}"
                )
        
                #Step5
                shell("cat {U3r1} {O4r1} > {OFr1}")
                shell("cat {U3r2} {O4r2} > {OFr2}")


elif config["libtype"]=="SA":
    rule cutadapt_SA:
        """Trims given paired-end reads with given parameters"""
        input:
            read1=expand("UMI_fastq/{sample}_R1.fq.gz", sample=config["samples"]),
            read2=expand("UMI_fastq/{sample}_R2.fq.gz", sample=config["samples"]),        
        output:
            finaloutread1=expand("cutadapt/{sample}_1.trimmed.fq.gz", sample=config["samples"]),
            finaloutread2=expand("cutadapt/{sample}_2.trimmed.fq.gz", sample=config["samples"]),
            statsouts=expand("cutadapt/{sample}.cutadaptstats.txt", sample=config["samples"]),
        threads: config["threads"]
        run:
            for Ir1,Ir2,S1,OFr1,OFr2 in zip(input.read1,input.read2,output.statsouts,output.finaloutread1,output.finaloutread2):
                shell(
                    "cutadapt -u 26 -U 26 --minimum-length 25 "
                    "-j {threads} -o {OFr1} -p {OFr2} {Ir1} {Ir2} > {S1}"
                )


else:
    print("libtype must be either 'LEXO' or 'SA'")

