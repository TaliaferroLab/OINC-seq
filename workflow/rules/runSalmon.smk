"""
Rules for transcriptome alignment with Salmon
For usage, include this in your workflow.
"""

if not config["dedupUMI"] and config["libtype"]=="LEXO":
    rule runSalmon:
        """aligns trimmed reads to transcriptome with Salmon"""
        input:
            read1=expand("PIGPEN_alignments/{sample}/STAR/{sample}.STARaligned.r1.fq.gz", sample=config["samples"]),
            read2=expand("PIGPEN_alignments/{sample}/STAR/{sample}.STARaligned.r2.fq.gz", sample=config["samples"]),
            Salmonindex=config["Salmonindex"],
        params:
            samplename=expand("PIGPEN_alignments/{sample}/salmon", sample=config["samples"]),
        output:
            sf=expand("PIGPEN_alignments/{sample}/salmon/{sample}.quant.sf", sample=config["samples"]),
            bam=temp(expand("PIGPEN_alignments/{sample}/salmon/{sample}.salmon.bam", sample=config["samples"])),
        threads: config["threads"]
        run:
            for r1,r2,samplename,bam in zip(input.read1,input.read2,params.samplename,output.bam):
                shell(
                    "salmon quant --libType A -p {threads} --seqBias --gcBias "
                    "--validateMappings -1 {r1} -2 {r2} -o {sample} "
                    "--index {input.Salmonindex} --writeMappings={bam} --writeQualities"
                )

elif config["dedupUMI"] and config["libtype"]=="LEXO":
    rule runSalmonDedup:
        """aligns trimmed reads to transcriptome with Salmon"""
        input:
            read1=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.r1.fq.gz", sample=config["samples"]),
            read2=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.r2.fq.gz", sample=config["samples"]),
            Salmonindex=config["Salmonindex"],
        params:
            samplename=expand("PIGPEN_alignments/{sample}/salmon", sample=config["samples"]),
        output:
            sf=expand("PIGPEN_alignments/{sample}/salmon/{sample}.quant.sf", sample=config["samples"]),
            bam=temp(expand("PIGPEN_alignments/{sample}/salmon/{sample}.salmon.bam", sample=config["samples"])),
        threads: config["threads"]
        run:
            for r1,r2,samplename,sf,bam in zip(input.read1,input.read2,params.samplename,output.sf,output.bam):
                shell(
                    "salmon quant --libType A -p {threads} --seqBias --gcBias "
                    "--validateMappings -1 {r1} -2 {r2} -o {samplename} "
                    "--index {input.Salmonindex} --writeMappings={bam} --writeQualities"
                )
                shell("mv {samplename}/quant.sf {sf}")
