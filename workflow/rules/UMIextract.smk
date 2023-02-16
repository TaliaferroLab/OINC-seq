"""
Rules for extraction of UMIs with UMI_tools
For usage, include this in your workflow.
"""

rule UMIextract:
    """extracts UMIs from raw fastq files"""
    input:
        read1=expand("RAWREADS/{sample}.R1.fq.gz", sample=config["samples"]),
        read2=expand("RAWREADS/{sample}.R2.fq.gz", sample=config["samples"]),
    output:
        read1=expand("UMI_fastq/{sample}_R1.fq.gz", sample=config["samples"]),
        read2=expand("UMI_fastq/{sample}_R2.fq.gz", sample=config["samples"]),
    run:
        if not config["libtype"]:
            print("libtype must be included in the config file")
        elif config["libtype"]:
            for r1,r2,o1,o2 in zip(input.read1,input.read2,output.read1,output.read2):
                if config["libtype"]=="LEXO":
                    shell(
                        "umi_tools extract -I {r1} "
                        "--bc-pattern=NNNNNN --read2-in={r2} "
                        "--stdout={o1} --read2-out={o2}"
                    )
                elif config["libtype"]=="SA":
                     shell(
                         "umi_tools extract -I {r2} "
                         "--bc-pattern=NNNNNNNNNNNN --read2-in={r1} "
                         "--stdout={o2} --read2-out={o1}"
                     )
                else:
                    print("libtype must be either 'LEXO' or 'SA'")
