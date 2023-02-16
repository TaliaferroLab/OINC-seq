"""
Rules for writing bam to paired fq files
For usage, include this in your workflow.
"""


if not config["dedupUMI"] and config["libtype"]=="LEXO":
    rule STARbam2FQ:
        """converts STAR aligned bam to paired fq"""
        input:
            bam=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam", sample=config["samples"]),
            bai=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam.bai", sample=config["samples"]),
        output:
            sortedbam=temp(expand("PIGPEN_alignments/{sample}/STAR/{sample}.temp.namesort.bam", sample=config["samples"])),
            fq1=temp(expand("PIGPEN_alignments/{sample}/STAR/{sample}.STARaligned.r1.fq.gz", sample=config["samples"])),
            fq2=temp(expand("PIGPEN_alignments/{sample}/STAR/{sample}.STARaligned.r2.fq.gz", sample=config["samples"])),
        threads: config["threads"]
        run:
            for bam,sortedbam,fq1,fq2 in zip(input.bam,output.sortedbam,output.fq1,output.fq2):
                shell(
                    "samtools collate --threads {threads} -u -o {sortedbam} {bam}"
                )
                shell(
                    "samtools fastq --threads {threads} -1 {fq1} -2 {fq2} "
                    "-0 /dev/null -s /dev/null -n {sortedbam}"
                )
    
elif config["dedupUMI"] and config["libtype"]=="LEXO":
    rule Dedupbam2FQ:
        """converts deduplicated STAR bam to paired fq"""
        input:
            bam=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam", sample=config["samples"]),
            bai=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam.bai", sample=config["samples"]),
        output:
            sortedbam=temp(expand("PIGPEN_alignments/{sample}/STAR/{sample}.temp.namesort.bam", sample=config["samples"])),
            fq1=temp(expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.r1.fq.gz", sample=config["samples"])),
            fq2=temp(expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.r2.fq.gz", sample=config["samples"])),
        threads: config["threads"]
        run:
            for bam,sortedbam,fq1,fq2 in zip(input.bam,output.sortedbam,output.fq1,output.fq2):
                shell(
                    "samtools collate --threads {threads} -u -o {sortedbam} {bam}"
                )
                shell(
                    "samtools fastq --threads {threads} -1 {fq1} -2 {fq2} "
                    "-0 /dev/null -s /dev/null -n {sortedbam}"
                )

