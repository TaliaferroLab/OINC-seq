"""
Rules for UMI deduplication with UMI-tools
For usage, include this in your workflow.
"""

if config["dedupUMI"]:
    rule dedupSTAR:
        """aligns trimmed reads to genome with STAR"""
        input:
            bam=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam", sample=config["samples"]),
            bai=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam.bai", sample=config["samples"]),

        output:
            bam=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam", sample=config["samples"]),
            bai=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam.bai", sample=config["samples"]),
        run:
            for Inbam,Outbam in zip(input.bam,output.bam):
                if config["libtype"]=="LEXO":
                    shell("umi_tools dedup -I {Inbam} --paired -S {Outbam}")
                    shell("samtools index {Outbam}")
                elif config["libtype"]=="SA":
                    shell("umi_tools dedup -I {Inbam} --paired --method=unique -S {Outbam}")
                    shell("samtools index {Outbam}")
                else:
                    print('LibType must be either "LEXO" or "SA".')

