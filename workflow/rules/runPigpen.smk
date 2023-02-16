"""
Rules for running pigpen.py
For usage, include this in your workflow.
"""

if config["libtype"]=="LEXO" and config["dedupUMI"]:
    rule runPigpen_LEXO_UMI:
        """runs Pigpen with desired parameters"""
        input:
            pigpen=config["pigpen"]["pigpen"],
            dedupbam=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam", sample=config["samples"]),
            dedupbai=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam.bai", sample=config["samples"]),
            PMbam=expand("PIGPEN_alignments/{sample}/postmaster/{sample}.postmaster.bam", sample=config["samples"]),
            PMbai=expand("PIGPEN_alignments/{sample}/postmaster/{sample}.postmaster.bam", sample=config["samples"]),
        params:
            samples=",".join(config["samples"]),
            gff=" --gff " + config["pigpen"]["gff"] if config["pigpen"]["gff"] else "",
            fa=" --genomeFasta " + config["pigpen"]["genomeFASTA"] if config["pigpen"]["genomeFASTA"] else "",
            CS=" --controlsamples " + config["pigpen"]["controlsamples"] if config["pigpen"]["controlsamples"] else "",
            VCF=" --snpfile " + config["pigpen"]["snpfile"] if config["pigpen"]["snpfile"] else "",
            MB=" --maskbed " + config["pigpen"]["maskbed"] if config["pigpen"]["maskbed"] else "",
            RB=" --ROIbed " + config["pigpen"]["ROIbed"] if config["pigpen"]["ROIbed"] else "",
            SC=" --SNPcoverage " + config["pigpen"]["SNPcoverage"] if config["pigpen"]["SNPcoverage"] else "",
            SF=" --SNPfreq " + config["pigpen"]["SNPfreq"] if config["pigpen"]["SNPfreq"] else "",
            NC=" --nConv " + config["pigpen"]["nconv"] if config["pigpen"]["nconv"] else "",
            MQ=" --minMappingQual " + config["pigpen"]["minMappingQual"] if config["pigpen"]["minMappingQual"] else "",
            OT=" " + config["pigpen"]["tags"] if config["pigpen"]["tags"] else "",
            OD=" --outputDir " + config["pigpen"]["outputDir"] if config["pigpen"]["outputDir"] else "",
        output:
            expand("PIGPEN_alignments/{outDir}/{sample}.pigpen.txt", sample=config["samples"], outDir =config["pigpen"]["outputDir"]),
        threads: config["threads"]
        run:
            os.chdir('PIGPEN_alignments')
            shell(
                "python -u {input.pigpen} --samplenames {params.samples} --nproc {threads}"
                "{params.gff}{params.fa}{params.CS}{params.VCF}{params.MB}{params.RB}"
                "{params.SC}{params.SF}{params.NC}{params.MQ}{params.OT}{params.OD}"
            )
            
elif config["libtype"]=="LEXO" and not config["dedupUMI"]:
    rule runPigpen_LEXO:
        """runs Pigpen with desired parameters"""
        input:
            pigpen=config["pigpen"]["pigpen"],
            STARbam=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam", sample=config["samples"]),
            STARbai=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam.bai", sample=config["samples"]),
            PMbam=expand("PIGPEN_alignments/{sample}/postmaster/{sample}.postmaster.bam", sample=config["samples"]),
            PMbai=expand("PIGPEN_alignments/{sample}/postmaster/{sample}.postmaster.bam", sample=config["samples"]),
        params:
            samples=",".join(config["samples"]),
            gff=" --gff " + config["pigpen"]["gff"] if config["pigpen"]["gff"] else "",
            fa=" --genomeFasta " + config["pigpen"]["genomeFASTA"] if config["pigpen"]["genomeFASTA"] else "",
            CS=" --controlsamples " + config["pigpen"]["controlsamples"] if config["pigpen"]["controlsamples"] else "",
            VCF=" --snpfile " + config["pigpen"]["snpfile"] if config["pigpen"]["snpfile"] else "",
            MB=" --maskbed " + config["pigpen"]["maskbed"] if config["pigpen"]["maskbed"] else "",
            RB=" --ROIbed " + config["pigpen"]["ROIbed"] if config["pigpen"]["ROIbed"] else "",
            SC=" --SNPcoverage " + config["pigpen"]["SNPcoverage"] if config["pigpen"]["SNPcoverage"] else "",
            SF=" --SNPfreq " + config["pigpen"]["SNPfreq"] if config["pigpen"]["SNPfreq"] else "",
            NC=" --nConv " + config["pigpen"]["nconv"] if config["pigpen"]["nconv"] else "",
            MQ=" --minMappingQual " + config["pigpen"]["minMappingQual"] if config["pigpen"]["minMappingQual"] else "",
            OT=" " + config["pigpen"]["tags"] if config["pigpen"]["tags"] else "",
            OD=" --outputDir " + config["pigpen"]["outputDir"] if config["pigpen"]["outputDir"] else "",
        output:
            expand("PIGPEN_alignments/{outDir}/{sample}.pigpen.txt", sample=config["samples"], outDir =config["pigpen"]["outputDir"]),
        threads: config["threads"]
        run:
            os.chdir('PIGPEN_alignments')
            shell(
                "python -u {input.pigpen} --samplenames {params.samples} --nproc {threads}"
                "{params.gff}{params.fa}{params.CS}{params.VCF}{params.MB}{params.RB}"
                "{params.SC}{params.SF}{params.NC}{params.MQ}{params.OT}{params.OD}"
            )

elif config["libtype"]=="SA" and config["dedupUMI"]:
    rule runPigpen_SA_UMI:
        """runs Pigpen with desired parameters"""
        input:
            pigpen=config["pigpen"]["pigpen"],
            dedupbam=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam", sample=config["samples"]),
            dedupbai=expand("PIGPEN_alignments/{sample}/STAR/{sample}.dedup.bam.bai", sample=config["samples"]),
        params:
            samples=",".join(config["samples"]),
            gff=" --gff " + config["pigpen"]["gff"] if config["pigpen"]["gff"] else "",
            fa=" --genomeFasta " + config["pigpen"]["genomeFASTA"] if config["pigpen"]["genomeFASTA"] else "",
            CS=" --controlsamples " + config["pigpen"]["controlsamples"] if config["pigpen"]["controlsamples"] else "",
            VCF=" --snpfile " + config["pigpen"]["snpfile"] if config["pigpen"]["snpfile"] else "",
            MB=" --maskbed " + config["pigpen"]["maskbed"] if config["pigpen"]["maskbed"] else "",
            RB=" --ROIbed " + config["pigpen"]["ROIbed"] if config["pigpen"]["ROIbed"] else "",
            SC=" --SNPcoverage " + config["pigpen"]["SNPcoverage"] if config["pigpen"]["SNPcoverage"] else "",
            SF=" --SNPfreq " + config["pigpen"]["SNPfreq"] if config["pigpen"]["SNPfreq"] else "",
            NC=" --nConv " + config["pigpen"]["nconv"] if config["pigpen"]["nconv"] else "",
            MQ=" --minMappingQual " + config["pigpen"]["minMappingQual"] if config["pigpen"]["minMappingQual"] else "",
            OT=" " + config["pigpen"]["tags"] if config["pigpen"]["tags"] else "",
            OD=" --outputDir " + config["pigpen"]["outputDir"] if config["pigpen"]["outputDir"] else "",
        output:
            expand("PIGPEN_alignments/{outDir}/{sample}.pigpen.txt", sample=config["samples"], outDir =config["pigpen"]["outputDir"]),
        threads: config["threads"]
        run:
            os.chdir('PIGPEN_alignments')
            shell(
                "python -u {input.pigpen} --samplenames {params.samples} --nproc {threads}"
                "{params.gff}{params.fa}{params.CS}{params.VCF}{params.MB}{params.RB}"
                "{params.SC}{params.SF}{params.NC}{params.MQ}{params.OT}{params.OD}"
            )

elif config["libtype"]=="SA" and not config["dedupUMI"]:
    rule runPigpen_SA:
        """runs Pigpen with desired parameters"""
        input:
            pigpen=config["pigpen"]["pigpen"],
            STARbam=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam", sample=config["samples"]),
            STARbai=expand("PIGPEN_alignments/{sample}/STAR/{sample}Aligned.sortedByCoord.out.bam.bai", sample=config["samples"]),
        params:
            samples=",".join(config["samples"]),
            gff=" --gff " + config["pigpen"]["gff"] if config["pigpen"]["gff"] else "",
            fa=" --genomeFasta " + config["pigpen"]["genomeFASTA"] if config["pigpen"]["genomeFASTA"] else "",
            CS=" --controlsamples " + config["pigpen"]["controlsamples"] if config["pigpen"]["controlsamples"] else "",
            VCF=" --snpfile " + config["pigpen"]["snpfile"] if config["pigpen"]["snpfile"] else "",
            MB=" --maskbed " + config["pigpen"]["maskbed"] if config["pigpen"]["maskbed"] else "",
            RB=" --ROIbed " + config["pigpen"]["ROIbed"] if config["pigpen"]["ROIbed"] else "",
            SC=" --SNPcoverage " + config["pigpen"]["SNPcoverage"] if config["pigpen"]["SNPcoverage"] else "",
            SF=" --SNPfreq " + config["pigpen"]["SNPfreq"] if config["pigpen"]["SNPfreq"] else "",
            NC=" --nConv " + config["pigpen"]["nconv"] if config["pigpen"]["nconv"] else "",
            MQ=" --minMappingQual " + config["pigpen"]["minMappingQual"] if config["pigpen"]["minMappingQual"] else "",
            OT=" " + config["pigpen"]["tags"] if config["pigpen"]["tags"] else "",
            OD=" --outputDir " + config["pigpen"]["outputDir"] if config["pigpen"]["outputDir"] else "",
        output:
            expand("PIGPEN_alignments/{outDir}/{sample}.pigpen.txt", sample=config["samples"], outDir =config["pigpen"]["outputDir"]),
        threads: config["threads"]
        run:
            os.chdir('PIGPEN_alignments')
            shell(
                "python -u {input.pigpen} --samplenames {params.samples} --nproc {threads}"
                "{params.gff}{params.fa}{params.CS}{params.VCF}{params.MB}{params.RB}"
                "{params.SC}{params.SF}{params.NC}{params.MQ}{params.OT}{params.OD}"
            )             
                

