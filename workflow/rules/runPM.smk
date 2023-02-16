"""
Rules for alignment proability extraction with postmaster
For usage, include this in your workflow.
"""

if config["libtype"]=="LEXO":
    rule runPostmaster:
        """extacts alignment probabilities from salmon bams and quantification"""
        input:
            sf=expand("PIGPEN_alignments/{sample}/salmon/{sample}.quant.sf", sample=config["samples"]),
            bam=expand("PIGPEN_alignments/{sample}/salmon/{sample}.salmon.bam", sample=config["samples"]),
        output:
            bam=expand("PIGPEN_alignments/{sample}/postmaster/{sample}.postmaster.bam", sample=config["samples"]),
            bai=expand("PIGPEN_alignments/{sample}/postmaster/{sample}.postmaster.bam.bai", sample=config["samples"]),
        threads: config["threads"]
        run:
            for sf,Sbam,PMbam in zip(input.sf,input.bam,output.bam):
                shell(
                    "postmaster --num-threads {threads} --quant {sf} "
                    "--alignments {Sbam} --output {PMbam}"
                )
                shell("samtools sort -@ {threads} -o {PMbam} {PMbam}")
                shell("samtools index {PMbam}")

