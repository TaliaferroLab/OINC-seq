"""
Rules for genome alignment with STAR
For usage, include this in your workflow.
"""

rule runSTAR:
    """aligns trimmed reads to genome with STAR"""
    input:
        read1=expand("cutadapt/{sample}_1.trimmed.fq.gz", sample=config["samples"]),
        read2=expand("cutadapt/{sample}_2.trimmed.fq.gz", sample=config["samples"]),
        STARindex=config["STARindex"],
    params:
        outprefix=expand("PIGPEN_alignments/{sample}/STAR/{sample}", sample=config["samples"]),
    output:
        bam=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam", sample=config["samples"]),
        bai=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam.bai", sample=config["samples"]),
    threads: config["threads"]
    run:
        for r1,r2,prefix,bam in zip(input.read1,input.read2,params.outprefix,output.bam):
            shell(
                "STAR --runMode alignReads --runThreadN {threads} "
                "--genomeLoad NoSharedMemory --genomeDir {input.STARindex} "
                "--readFilesIn {r1} {r2} --readFilesCommand zcat "
                "--outFileNamePrefix {prefix} --outSAMtype BAM "
                "SortedByCoordinate --outSAMstrandField intronMotif "
                "--outSAMattributes MD NH --outSAMmultNmax 1"
            )
            shell("samtools index {bam}")
            
